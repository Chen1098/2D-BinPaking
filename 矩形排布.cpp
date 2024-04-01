#include <iostream>
#include <vector>
#include<bits/stdc++.h>
using namespace std;
struct Rectangle {
	float width, height, x, y;
	int id;
};

struct RGB{
	int r;
	int g;
	int b;
} rgb;

RGB color(){
	RGB use;
	use.r=(rand()%255)*rand()%100%255;
	use.g=rand()%255;
	use.b=rand()%255;
	return use;
}

bool doesOverlap(const Rectangle& a, const Rectangle& b) {
	return (a.x < b.x + b.width && a.x + a.width > b.x &&
	        a.y < b.y + b.height && a.y + a.height > b.y);
}

bool isWithinBoard(const Rectangle& r, float boardWidth, float boardHeight) {
	return r.x + r.width <= boardWidth && r.y + r.height <= boardHeight;
}

void placeRectangles(std::vector<Rectangle>& rectangles, float boardWidth, float boardHeight) {
	for(size_t i = 0; i < rectangles.size(); i++) {
		bool placed = false;
		while (!placed && isWithinBoard(rectangles[i], boardWidth, boardHeight)) {
			placed = true; // Assume we've placed until proven otherwise.
			for (size_t j = 0; j < i; j++) {
				if (doesOverlap(rectangles[i], rectangles[j])) {
					rectangles[i].x += 1.0; // move it to the right
					if (rectangles[i].x + rectangles[i].width > boardWidth) { // if it goes out of board width, reset x and shift it down
						rectangles[i].x = 0;
						rectangles[i].y += 1.0;
					}
					placed = false; // we've moved the rectangle, so we need to check placement again
					break;
				}
			}
		}
		if (!isWithinBoard(rectangles[i], boardWidth, boardHeight)) {
			std::cout << "矩形 " << i+1 << "放不下" << std::endl;
		}
	}
}
struct DOT{
	double x[1000],y[1000];
	int cnt;
} shapes[10001];
Rectangle change_rec(int nowat){
	double maxx=-0x3f3f3f3f,maxy=-0x3f3f3f3f,minx=0x3f3f3f3f,miny=0x3f3f3f3f;
	for(int i=1;i<=shapes[nowat].cnt;i++){
		maxx=fmax(maxx,shapes[nowat].x[i]);
		maxy=fmax(maxy,shapes[nowat].y[i]);
		minx=fmin(minx,shapes[nowat].x[i]);
		miny=fmin(miny,shapes[nowat].y[i]);
	}
	Rectangle now;
	now.height=maxy-miny;
	now.width=maxx-minx;
	now.x=0;
	now.y=0;
	now.id=nowat;
	return now;
}
int main() {
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout); 
	int in_cnt;
	float boardWidth;
	float boardHeight;
	std::vector<Rectangle> rectangles;
	cout<<"输入底板的长与宽"<<endl;
	cin>>boardWidth>>boardHeight; 
	cout << "输入图形的个数 "<<endl;
	cin>>in_cnt;
	for(int i=1;i<=in_cnt;i++){
		Rectangle use;
		int sides;
		cout<<"请输入图形边的个数"<<endl;
		cin>>shapes[i].cnt;
		for(int j=1;j<=shapes[i].cnt;j++){
			cin>>shapes[i].x[j]>>shapes[i].y[j];
		} 
		use=change_rec(i);
		rectangles.push_back(use);
	}


	placeRectangles(rectangles, boardWidth, boardHeight);

	for(const auto& rect : rectangles) {
		if (isWithinBoard(rect, boardWidth, boardHeight)) {
			std::cout << "Rectangle at (" << rect.x << "," << rect.y << ")" << std::endl;
			
		}
	}
	cout<<"输出矩形"<<endl<<endl; 
	cout<<"<!DOCTYPE html>"<<endl<<"<html>"<<endl<<"<body>"<<endl;
	cout<<"<svg  height=\"9999\" width=\"9999\">"<<endl;
	for (const auto& rect : rectangles) {
		cout<<"<polygon points="<<"\"";
		cout << rect.x*50 << "," << rect.y*50 <<' '<<rect.x*50+rect.width*50<<','<<rect.y*50<<' '<<rect.x*50+rect.width*50<<','<<rect.y*50+rect.height*50<<' '<<rect.x*50<<','<<rect.y*50+rect.height*50<<endl;
		cout<<"\""<<"style=\""<<"fill:rgb("<<color().r<<','<<color().b<<','<<color().g<<");stroke:purple;stroke-width:1"<<"\""<<"/>"<<endl;
	}
	cout<<"</svg>"<<"</body>"<<"</html>"<<endl;
	cout<<"输出不规则图形"<<endl<<endl;
	
	
	cout<<"<!DOCTYPE html>"<<endl<<"<html>"<<endl<<"<body>"<<endl;
	cout<<"<svg  height=\"9999\" width=\"9999\">"<<endl;
	for(const auto& rect:rectangles){
		cout<<"<polygon points="<<"\"";
		for(int i=1;i<=shapes[rect.id].cnt;i++){
			cout<<rect.x*50+shapes[rect.id].x[i]*50<<","<<rect.y*50+shapes[rect.id].y[i]*50<<" ";
		}
		cout<<"\""<<"style=\""<<"fill:rgb("<<color().r<<','<<color().b<<','<<color().g<<");stroke:purple;stroke-width:1"<<"\""<<"/>"<<endl;
	} 
	return 0;
}
/*
20 20
9
7 3
3 5
3 5
3 7
3 2
3 7
3 2
3 1
8 3
*/

/*
30 30
23
4
0 0
0 7
3 7
3 0
4
0 0
0 3
5 3
5 0
4
0 0
0 3
5 3
5 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
2 3
2 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
2 3
2 0
4
0 0
0 3
1 3
1 0
4
0 0
0 3
8 3
8 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
5 3
5 0
4
0 0
0 3
5 3
5 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
2 3
2 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
2 3
2 0
4
0 0
0 3
1 3
1 0
4
0 0
0 3
8 3
8 0
4
0 0
0 7
3 7
3 0
4
0 0
0 3
2 3
2 0
4
0 0
0 3
1 3
1 0
4
0 0
0 3
8 3
8 0
3
0 0
2 0
2 2
*/
