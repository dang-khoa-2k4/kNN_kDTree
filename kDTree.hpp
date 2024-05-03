#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */

extern size_t ELEMENT_INDEX;

struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    int label;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }
    kDTreeNode(vector<int> data, int label, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
        this->label = label;
    }
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    int size;
    struct ListWithLabel
    {
        vector<int> data;
        int label;
        ListWithLabel(vector<int> data, int label)
        {
            this->data = data;
            this->label = label;
        }
        bool operator<(const ListWithLabel &other) const
        {
            return data < other.data;
        }
    };

    // some drapper function 
    int height(kDTreeNode * root) const;
    double distance(const vector<int> &v1, const vector<int> &v2) const;
    int leafCount(kDTreeNode * root) const;
    void traverse(kDTreeNode * root, void (*f)(kDTreeNode *)) const;
    void traverse(kDTreeNode * root, kDTreeNode * newNode, void (*f)(kDTreeNode * curNode, kDTreeNode * newNode)) const;
    void inorderTraversal(kDTreeNode * root) const;
    void preorderTraversal(kDTreeNode * root) const;
    void postorderTraversal(kDTreeNode * root) const;
    void remove(kDTreeNode * &root, const vector<int> &point, int axis);
    void buildTree(kDTreeNode * &root, vector<vector<int>> &pointList, int axis);
    void buildTree(kDTreeNode * &root, vector<ListWithLabel> &pointList, int axis);
    void nearestNeighbour(kDTreeNode * root, const vector<int> &target, kDTreeNode *best, int axis);
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
    void insert(const vector<int> &point, int label);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void buildTree(const vector<vector<int>> &pointList, vector<int> &labelList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
    // extra function 
    // DFS preorder
    void traverse(void (*f)(kDTreeNode *)) const;  
    void traverse(kDTreeNode * newNode, void (*f)(kDTreeNode *, kDTreeNode *)) const;  
    kDTreeNode *findMin(kDTreeNode * root, int d, int axis);
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    kDTree *treeData;
    int numClasses;
    std::vector<std::vector<int>> convertListListToVectorVector(const std::list<std::list<int>>& list_of_lists);
    template<typename T> std::list<T> convertVectorToList(const std::vector<T>& input_list);
    template<typename T> std::vector<T> convertListToVector(const std::list<T>& input_list);
public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
inline bool operator<(const std::vector<int> &v1, const std::vector<int> &v2) {
    return v1[ELEMENT_INDEX] < v2[ELEMENT_INDEX];
}
template<class T> void merge(T * arr, int l, int m, int r);
template<class T> void mergeSort(T * arr, int l, int r);