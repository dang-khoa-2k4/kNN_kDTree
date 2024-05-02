#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    int size;

    // some drapper function 
    int height(kDTreeNode * root) const;
    int leafCount(kDTreeNode * root) const;
    void traverse(kDTreeNode * root, void (*f)(kDTreeNode *)) const;
    void traverse(kDTreeNode * root, kDTreeNode * newNode, void (*f)(kDTreeNode * curNode, kDTreeNode * newNode)) const;
    void inorderTraversal(kDTreeNode * root) const;
    void preorderTraversal(kDTreeNode * root) const;
    void postorderTraversal(kDTreeNode * root) const;
public:
    kDTree(int k = 2);
    ~kDTree();

    const kDTree &operator=(const kDTree &other);
    kDTree(const kDTree &other);

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const;
    int nodeCount() const;
    int leafCount() const;
    
    void insert(const vector<int> &point);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
    // extra function 
    // DFS postorder
    void traverse(void (*f)(kDTreeNode *)) const;  
    void traverse(kDTreeNode * newNode, void (*f)(kDTreeNode *, kDTreeNode *)) const;  
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    int numClasses;

public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed