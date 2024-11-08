#include "rayAccelerator.h"
#include "macros.h"

using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
			//this->n_objs = n_objs_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }


void BVH::Build(vector<Object *> &objs) {

		
			BVHNode *root = new BVHNode();

			Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			AABB world_bbox = AABB(min, max);

			for (Object* obj : objs) {
				AABB bbox = obj->GetBoundingBox();
				world_bbox.extend(bbox);
				objects.push_back(obj);
			}
			world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
			world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
			root->setAABB(world_bbox);
			nodes.push_back(root);
			build_recursive(0, objects.size(), root); // -> root node takes all the 
		}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {

		//right_index, left_index and split_index refer to the indices in the objects vector
	   // do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	    // node.index can have a index of objects vector or a index of nodes vector
			
		// Calculate the AABB for the objects in the specified range
		AABB bounding_box;
		for (int i = left_index; i < right_index; i++) {
			Object* obj = objects[i];
			AABB bbox = obj->GetBoundingBox();
			bounding_box.extend(bbox); //emcompasses all objects
		}
		node->setAABB(bounding_box); //sets bb to current node
	
		// Check if the number of objects in the range is below a threshold or if the recursion has reached a maximum depth
		int num_objs = right_index - left_index;
		// recursively analyzes each node until each node contains two objects most likely to collide.
		if (num_objs <= Threshold) {
			// Create a leaf node and store the range of objects
			node->makeLeaf(left_index, num_objs);
			return;
		}

		// Calculate the split position = mid point
		int split_index = (left_index + right_index) / 2;

		//longest axis
		AABB node_bb = node->getAABB();

		int split_axis;

		Vector dist = node_bb.max - node_bb.min;

		if (dist.x >= dist.y && dist.x >= dist.z) {
			split_axis = 0; // X axis is the longest
		}
		else if (dist.y >= dist.x && dist.y >= dist.z) {
			split_axis = 1; // Y axis
		}
		else {
			split_axis = 2; // Z axis
		}

		std::sort(objects.begin() + left_index, objects.begin() + right_index,
			[split_axis](Object* a, Object* b) {
				Vector centroid_a = a->getCentroid();
				Vector centroid_b = b->getCentroid();

				return centroid_a.getAxisValue(split_axis) < centroid_b.getAxisValue(split_axis);
			});

		// Create a node and recursively build the left and right child nodes
		BVHNode* left_node = new BVHNode();
		BVHNode* right_node = new BVHNode();

		nodes.push_back(left_node);
		nodes.push_back(right_node);

		node->makeNode(nodes.size() - 2);

		build_recursive(left_index, split_index, left_node);
		build_recursive(split_index, right_index, right_node);
	}

// find the closest intersection of the ray with the objects in the BVH tree.
bool BVH::Traverse(Ray& ray, Object** hit_obj, Vector& hit_point) {
			float tmp;
			float tmin = FLT_MAX;  //contains the closest primitive intersection
			bool hit = false;

			BVHNode* currentNode = nodes[0];

			// Traversal stack for storing intermediate nodes
			std::stack<BVHNode*> traversalStack;
			traversalStack.push(currentNode);

			//iterates through the nodes in the BVH tree, checking for intersection between the ray and the current node's bounding box
			while (!traversalStack.empty()) {
				currentNode = traversalStack.top();
				traversalStack.pop();

				// If the ray does not intersect with the current bounding box, skip this node
				if (!currentNode->getAABB().intercepts(ray,tmp))
					continue;

				// If the node is a leaf node, check for intersections between the ray and the objects contained in the leaf
				if (currentNode->isLeaf()) {
					// Leaf node, iterate through the objects in the leaf and find the closest intersection
					for (unsigned int i = currentNode->getIndex(); i < currentNode->getIndex() + currentNode->getNObjs(); ++i) {
						// If an intersection is found and it is closer than the previous closest
						Object* currentObj = objects[i];

						//printf("intercepts");
						if (currentObj->intercepts(ray, tmp) && tmp < tmin) {
							// Update the closest intersection
							
							tmin = tmp;
							*hit_obj = currentObj;
							//printf("hit point: %f, %f, %f\n", currentObj->getCentroid().x, currentObj->getCentroid().y, currentObj->getCentroid().z);
							
							hit_point = ray.origin + ray.direction * tmp;
							hit = true;
						}
					}
				}
				else {
					// Non-leaf node, traverse down the tree
					traversalStack.push(nodes[currentNode->getIndex()]);  // Push the right child node
					traversalStack.push(nodes[currentNode->getIndex() + 1]);  // Push the left child node
				}
			}

			return hit;
}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
			float tmp;

			double length = ray.direction.length(); //distance between light and intersection point
			ray.direction.normalize();

			BVHNode* currentNode = nodes[0];

			// Traversal stack for storing intermediate nodes
			std::stack<BVHNode*> traversalStack;
			traversalStack.push(currentNode);

			while (!traversalStack.empty()) {
				currentNode = traversalStack.top();
				traversalStack.pop();

				// If the ray does not intersect with the current bounding box, skip this node
				if (!currentNode->getAABB().intercepts(ray,tmp))
					continue;

				if (currentNode->isLeaf()) {
					// Leaf node, check for any intersection
					for (unsigned int i = currentNode->getIndex(); i < currentNode->getIndex() + currentNode->getNObjs(); ++i) {
						float t;
						if (objects[i]->intercepts(ray, t) && t < length)
							return true;  // Shadow ray intersects an object
					}
				}
				else {
					// Non-leaf node, traverse down the tree
					traversalStack.push(nodes[currentNode->getIndex()]);  // Push the right child node
					traversalStack.push(nodes[currentNode->getIndex() + 1]);  // Push the left child node
				}
			}

			return false; // Shadow ray does not intersect any object

}		
