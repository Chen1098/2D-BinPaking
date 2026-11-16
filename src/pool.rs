//! A tiny fork-join pool for very short parallel sections.
//!
//! One separation iteration is only a few hundred microseconds of work per worker, so the
//! wake-up latency of a general purpose pool costs a large fraction of the wall time.  Here
//! the helper threads spin briefly between jobs (and park when the pause gets long), and the
//! calling thread takes part as participant 0.

use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::Duration;

type Job = *const (dyn Fn(usize) + Sync);

struct Shared {
    epoch: AtomicU64,
    done: AtomicUsize,
    stop: AtomicBool,
    job: UnsafeCell<Option<Job>>,
}

unsafe impl Sync for Shared {}
unsafe impl Send for Shared {}

pub struct SpinPool {
    shared: Arc<Shared>,
    handles: Vec<JoinHandle<()>>,
    n: usize,
}

impl SpinPool {
    /// A pool with `n` participants (the caller plus `n - 1` helper threads).
    pub fn new(n: usize) -> Self {
        let n = n.max(1);
        let shared = Arc::new(Shared { epoch: AtomicU64::new(0), done: AtomicUsize::new(0), stop: AtomicBool::new(false), job: UnsafeCell::new(None) });
        let mut handles = Vec::with_capacity(n - 1);
        for idx in 1..n {
            let sh = shared.clone();
            handles.push(std::thread::spawn(move || {
                let mut seen = 0u64;
                loop {
                    let mut spins = 0u32;
                    loop {
                        let e = sh.epoch.load(Ordering::Acquire);
                        if e != seen {
                            seen = e;
                            break;
                        }
                        spins += 1;
                        if spins < 4000 {
                            std::hint::spin_loop();
                        } else if spins < 4200 {
                            std::thread::yield_now();
                        } else {
                            std::thread::park_timeout(Duration::from_millis(2));
                        }
                    }
                    if sh.stop.load(Ordering::Acquire) {
                        return;
                    }
                    let job = unsafe { (*sh.job.get()).expect("job set before the epoch changes") };
                    unsafe { (*job)(idx) };
                    sh.done.fetch_add(1, Ordering::Release);
                }
            }));
        }
        SpinPool { shared, handles, n }
    }

    pub fn size(&self) -> usize {
        self.n
    }

    /// Run `f(k)` for every participant index `k` in `0..n` and wait for all of them.
    pub fn run<F: Fn(usize) + Sync>(&self, f: F) {
        if self.n == 1 {
            f(0);
            return;
        }
        let fref: &(dyn Fn(usize) + Sync) = &f;
        // SAFETY: the pointer is only dereferenced by the helpers before `done` reaches
        // n - 1, and this function does not return (so `f` stays alive) until then.
        let ptr: Job = unsafe { std::mem::transmute::<&(dyn Fn(usize) + Sync), Job>(fref) };
        unsafe {
            *self.shared.job.get() = Some(ptr);
        }
        self.shared.done.store(0, Ordering::Release);
        self.shared.epoch.fetch_add(1, Ordering::Release);
        for h in &self.handles {
            h.thread().unpark();
        }
        f(0);
        let mut spins = 0u32;
        while self.shared.done.load(Ordering::Acquire) < self.n - 1 {
            spins += 1;
            if spins < 20000 {
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
        unsafe {
            *self.shared.job.get() = None;
        }
    }
}

impl Drop for SpinPool {
    fn drop(&mut self) {
        self.shared.stop.store(true, Ordering::Release);
        self.shared.epoch.fetch_add(1, Ordering::Release);
        for h in &self.handles {
            h.thread().unpark();
        }
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}
