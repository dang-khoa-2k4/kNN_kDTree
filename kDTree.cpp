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
    int axis = 0;
    kDTreeNode **pT = &root;
    while (*pT)
    {
        if (((*pT)->data)[axis % k] > point[axis % k]) pT = &((*pT)->left);
        else pT = &((*pT)->right);
        axis++;
    }
    *pT = new kDTreeNode(point);
    size++;
}

kDTreeNode *kDTree::findMin(kDTreeNode * root, int d, int axis)
{
    if (!root) return nullptr;
    if (axis % k == d)
    {
        if (!root->left) return root;
        return findMin(root->left, d, axis + 1);
    }
    kDTreeNode *left = findMin(root->left, d, axis + 1);
    kDTreeNode *right = findMin(root->right, d, axis + 1);
    kDTreeNode *minNode = (left->data[d] < root->data[d]) ? left : root;
    minNode = (right->data[d] < minNode->data[d]) ? right : minNode;
    return minNode;
}


void kDTree::remove(kDTreeNode * &root, const vector<int> &point, int axis)
{
    if (!root) return;
    if (root->data == point)
    {
        if (root->right)
        {
            kDTreeNode *pTemp = findMin(root->right, axis, axis + 1);
            root->data = pTemp->data;
            remove(root->right, pTemp->data, axis + 1);
            return;
        }
        else if (root->left)
        {
            kDTreeNode *pTemp = findMin(root->left, axis, axis + 1);
            root->data = pTemp->data;
            root->right = root->left;
            root->left = nullptr;
            remove(root->right, pTemp->data, axis + 1);
            return;
        }
        else
        {
            delete root;
            root = nullptr;
            size--;
        } 
    }
    else if ((root->data)[axis % k] > point[axis % k])
        remove(root->left, point, axis + 1);
    else
        remove(root->right, point, axis + 1);
}

void kDTree::remove(const vector<int> &point)
{
    return remove(root, point, 0);
}

bool kDTree::search(const vector<int> &point)
{
    int axis = 0;
    kDTreeNode *pT = root;
    while (pT)
    {
        if (pT->data == point) return true;
        if ((pT->data)[axis % k] > point[axis % k]) pT = pT->left;
        else pT = pT->right;
        axis++;
    }
    return false;
}

void kDTree::buildTree(kDTreeNode *&root, 
                       const std::vector<std::vector<int>> &pointList, int axis)
{
    if (pointList.empty()) return;
    if (pointList.size() == 1)
    {
        root = new kDTreeNode(pointList[0]);
        return;
    }
    int n = pointList.size();
    ELEMENT_INDEX = axis % k;
    mergeSort(&pointList, 0, n - 1);
    root = new kDTreeNode(pointList[n / 2]);
    std::vector<std::vector<int>> leftList(pointList.begin(), pointList.begin() + n / 2);
    std::vector<std::vector<int>> rightList(pointList.begin() + n / 2 + 1, pointList.end());
    buildTree(root->left, leftList, axis + 1);
    buildTree(root->right, rightList, axis + 1);
}

void kDTree::buildTree(const std::vector<std::vector<int>> &pointList)
{
    return buildTree(root, pointList, 0);
}

double kDTree::distance(const std::vector<int> &a, const std::vector<int> &b) const
{
    double sum = 0;
    for (int i = 0; i < (int) a.size(); i++)
    {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(sum);
}

void kDTree::nearestNeighbour(kDTreeNode *root, const std::vector<int> &target, kDTreeNode *best, int axis)
{
    if (!root) return;
    if ((root->data)[axis % k] < target[axis % k])
    {
        nearestNeighbour(root->right, target, best, axis + 1);
        if (!best || distance(root->data, target) < distance(best->data, target)) best = root;
        if (best && abs((root->data)[axis % k] - target[axis % k]) < distance(best->data, target))
            nearestNeighbour(root->left, target, best, axis + 1);
    }
    else
    {
        nearestNeighbour(root->left, target, best, axis + 1);
        if (!best || distance(root->data, target) < distance(best->data, target)) best = root;
        if (best && abs((root->data)[axis % k] - target[axis % k]) < distance(best->data, target))
            nearestNeighbour(root->right, target, best, axis + 1);
    }
}

void kDTree::nearestNeighbour(const std::vector<int> &target, kDTreeNode *best)
{
    return nearestNeighbour(root, target, best, 0);
}
void kDTree::kNearestNeighbour(const std::vector<int> &target, int k, std::vector<kDTreeNode *> &bestList)
{
    while (k--)
    {
        kDTreeNode *best = nullptr;
        nearestNeighbour(target, best);
        bestList.push_back(best);
        remove(best->data);
    }
}

// kNN part

kNN::kNN(int k)
{
    this->k = k;
    X_train = nullptr;
    y_train = nullptr;

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

template<class T>
void merge(T *arr, int l, int m, int r)
{
    T * tempL = new T[m - l + 1];
    for (int i = l; i < m + 1; i++)
    {
        tempL[i - l] = arr[i];
    }
    int i = 0, j = m + 1, k = l;
    while (i < m - l + 1 && j < r + 1)
    {
        if (tempL[i] < arr[j])
        {
            arr[k] = tempL[i];
            i++;
        }
        else 
        {
            arr[k] = arr[j];
            j++;
        }
        k++;
    }
    while (i < m - l + 1)
    {
        arr[k] = tempL[i];
        i++;
        k++;
    }
    delete[] tempL;
}

template<class T>
void mergeSort(T *arr, int l, int r)
{
    if (l < r)
    {
        int m = (l + r) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}