#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

// kDTree data structure part

kDTree::kDTree(int k)
{
    this->k = k;
    root = nullptr;
    size = 0;
}

kDTree::~kDTree()
{
    traverse([](kDTreeNode * nT) {  delete nT;  });
    k = 2;
    size = 0;
}

void kDTree::traverse(kDTreeNode * root, void (*f)(kDTreeNode *)) const
{
    if (!root) return;
    traverse(root->left, f);
    traverse(root->right, f);
    f(root);
}

void kDTree::traverse(kDTreeNode * root, kDTreeNode * newNode, void (*f)(kDTreeNode * curNode, kDTreeNode * newNode)) const
{
    if (!root) return;
    traverse(root->left, newNode->left, f);
    traverse(root->right, newNode->right, f);
    f(root, newNode);
}

void kDTree::traverse(void (*f)(kDTreeNode *)) const
{
    traverse(root, f);
}

void kDTree::traverse(kDTreeNode * newNode, void (*f)(kDTreeNode *, kDTreeNode *)) const
{
    traverse(root, newNode, f);
}

// traverse and create new tree
const kDTree &kDTree::operator=(const kDTree &other)
{
    kDTree * newTree = new kDTree(other.k);
    other.traverse(newTree->root,
    [](kDTreeNode * root, kDTreeNode *newNode) 
    {
        newNode = new kDTreeNode(root->data);
    });
    this->size = other.size;
    return *newTree;
}
// traverse and create new tree
kDTree::kDTree(const kDTree &other)
{
    other.traverse(root,
    [](kDTreeNode * root, kDTreeNode *newNode) 
    {
        newNode = new kDTreeNode(root->data);
    });
    this->k = other.k;
    this->size = other.size;
}   

void kDTree::inorderTraversal(kDTreeNode * root) const
{
    if (!root) return;
    inorderTraversal(root->left);
    cout << (root != this->root) ? " " : "";
    cout << '(' << (root->data)[0];
    for (int i = 1; i < (int) (root->data).size(); i++) 
        cout << ", " << (root->data)[i]; 
    cout << ')';
    inorderTraversal(root->right);
}
void kDTree::inorderTraversal() const
{
    if (!root) return;
    inorderTraversal(root);
}

void kDTree::preorderTraversal(kDTreeNode * root) const
{
    if (!root) return;
    cout << (root != this->root) ? " " : "";
    cout << '(' << (root->data)[0];
    for (int i = 1; i < (int) (root->data).size(); i++) 
        cout << ", " << (root->data)[i]; 
    cout << ')';
    preorderTraversal(root->left);
    preorderTraversal(root->right);
}
void kDTree::preorderTraversal() const
{
    if (!root) return;
    preorderTraversal(root);
}

void kDTree::postorderTraversal(kDTreeNode * root) const
{
    if (!root) return;
    postorderTraversal(root->left);
    postorderTraversal(root->right);
    cout << (root != this->root) ? " " : "";
    cout << '(' << (root->data)[0];
    for (int i = 1; i < (int) (root->data).size(); i++) 
        cout << ", " << (root->data)[i]; 
    cout << ')';
}
void kDTree::postorderTraversal() const
{
    if (!root) return;
    postorderTraversal(root);
}

int kDTree::height(kDTreeNode * root) const
{
    if (!root) return 0;
    int heightLeft = height(root->left);
    int heightRight = height(root->right);
    return 1 + (heightLeft > heightRight ? heightLeft : heightRight);
}

int kDTree::height() const 
{
    return height(root);
}

int kDTree::nodeCount() const
{
    return size;
}

int kDTree::leafCount(kDTreeNode * root) const
{
    if (!root) return 0;
    if (!root->left && !root->right) return 1;
    return leafCount(root->left) + leafCount(root->right);
}

int kDTree::leafCount() const
{
    return leafCount(root);
}

void kDTree::insert(const vector<int> &point)
{

}
void kDTree::remove(const vector<int> &point)
{

}
bool kDTree::search(const vector<int> &point)
{

}
void kDTree::buildTree(const std::vector<std::vector<int>> &pointList)
{

}
void kDTree::nearestNeighbour(const std::vector<int> &target, kDTreeNode *best)
{

}
void kDTree::kNearestNeighbour(const std::vector<int> &target, int k, std::vector<kDTreeNode *> &bestList)
{

}

// kNN part

kNN::kNN(int k)
{

}
void kNN::fit(Dataset &X_train, Dataset &y_train)
{

}
Dataset kNN::predict(Dataset &X_test)
{

}
double kNN::score(const Dataset &y_test, const Dataset &y_pred)
{

}
