#include "Node.h"

Node::Node()
{
	x = 0;
	y = 0;
	bc = 0;
}

Node::Node(float x, float y, bool bc)
{
	this->x = x;
	this->y = y;
	this->bc = bc;
}

void Node::set_x(float x)
{
	this->x = x;
}

void Node::set_y(float y)
{
	this->y = y;
}

void Node::set_bc(bool bc)
{
	this->bc = bc;
}

float Node::get_x()
{
	return x;
}

float Node::get_y()
{
	return y;
}

bool Node::get_bc()
{
	return bc;
}
