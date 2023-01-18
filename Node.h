#pragma once
class Node
{
	float x;
	float y;
	bool bc;
public:
	Node();
	Node(float x, float y, bool bc);
	void set_x(float x);
	void set_y(float y);
	void set_bc(bool bc);
	float get_x();
	float get_y();
	bool get_bc();
};

