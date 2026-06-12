// B-Tree
// description:

#pragma once

#include <vector>

// Note: btree also works for substrings, meaning that "2, 3" will be said to be
// in the tree if "1, 2, 3" is (it stores nodes in the reverse order of the
// vector) Implement btree for arbitrary containers. Also consider if using a
// list for children_values would be a better choice.
template <typename T> struct btree {
private:
  std::vector<T> childrenValues = {};
  std::vector<btree<T> *> childrenPointers = {};
  void updateNodeRecursively(std::vector<T> inputVector, int vectorSize) {
    if (vectorSize != 0) {
      vectorSize--;
      for (int i = 0; i < childrenValues.size(); i++) {
        if (childrenValues[i] == inputVector[vectorSize]) {
          (*childrenPointers[i]).updateNodeRecursively(inputVector, vectorSize);
          return;
        }
      }
      btree<T> *childTreePointer =
          new btree<T>; // need to use custom resource management so children
                        // don't get deleted when leaving "update" member
                        // function's scope; ADD CUSTOM DESTRUCTOR, AS WELL!
                        // CURRENTLY, ONLY ROOT GETS DELETED WHEN SCOPE ITS
                        // INITIATED IN ENDS
      childrenValues.push_back(inputVector[vectorSize]);
      childrenPointers.push_back(childTreePointer);
      (*childTreePointer).updateNodeRecursively(inputVector, vectorSize);
    }
  }
  bool searchNodeRecursively(std::vector<T> inputVector, int vectorSize) {
    if (vectorSize != 0) {
      vectorSize--;
      for (int i = 0; i < childrenValues.size(); i++) {
        if (childrenValues[i] == inputVector[vectorSize]) {
          return (*childrenPointers[i])
              .searchNodeRecursively(inputVector, vectorSize);
        }
      }
      return false;
    }
    return true;
  }

public:
  void insertVector(std::vector<T> inputVector) {
    updateNodeRecursively(inputVector, inputVector.size());
  }
  bool containsVector(std::vector<T> inputVector) {
    return searchNodeRecursively(inputVector, inputVector.size());
  }
  std::vector<T> getChildrenValues() { return childrenValues; }
  std::vector<btree<T> *> getChildrenPointers() { return childrenPointers; }
};