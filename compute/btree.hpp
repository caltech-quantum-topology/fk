// B-Tree
// description:

#pragma once

#include <vector>

// Note: btree also works for substrings, meaning that "2, 3" will be said to be in the tree if "1, 2, 3" is (it stores nodes in the reverse order of the vector)
// Implement btree for arbitrary containers. Also consider if using a list for children_values would be a better choice.
template <typename T>
struct btree {
    private:
        std::vector<T> children_values = {};
        std::vector<btree<T>*> children_pointers = {};
        void private_update (std::vector<T> v, int size) {
            if (size != 0) {
                size--;
                for (int i = 0; i < children_values.size(); i++) {
                    if (children_values[i] == v[size]) {
                        (*children_pointers[i]).private_update(v, size);
                        return;
                    }
                }
                btree<T>* btr_ptr = new btree<T>; // need to use custom resource management so children don't get deleted when leaving "update" member function's scope; ADD CUSTOM DESTRUCTOR, AS WELL! CURRENTLY, ONLY ROOT GETS DELETED WHEN SCOPE ITS INITIATED IN ENDS
                children_values.push_back(v[size]);
                children_pointers.push_back(btr_ptr);
                (*btr_ptr).private_update(v, size);
            }
        }
        bool private_contains (std::vector<T> v, int size) {
            if (size != 0) {
                size--;
                for (int i = 0; i < children_values.size(); i++) {
                    if (children_values[i] == v[size]) {
                        return (*children_pointers[i]).private_contains(v, size);
                    }
                }
                return false;
            }
            return true;
        }
    public:
        void update (std::vector<T> v) {
            private_update(v, v.size());
        }
        bool contains (std::vector<T> v) {
            return private_contains(v, v.size());
        }
        std::vector<T> get_children_values () {
            return children_values;
        }
        std::vector<btree<T>*> get_children_pointers () {
            return children_pointers;
        }
};