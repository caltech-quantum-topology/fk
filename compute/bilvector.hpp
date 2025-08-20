// "bi list vectors"
// description:

#pragma once

#include <vector>
#include <list>

template <typename T>
struct bilvector {
    private:
        int component_size;
        int nnvectors = 0;
        int npvectors = 0;
        int max_nindex = 0;
        int max_pindex = 0;
        T default_;
        std::list<std::vector<T>> nvectors = {};
        std::list<std::vector<T>> pvectors = {};
    public:
        bilvector (int initial_nnvectors, int initial_npvectors, int component_size_, T default__) {
            component_size = component_size_;
            default_ = default__;
            nnvectors = initial_nnvectors;
            for (int i = 0; i < initial_nnvectors; i++) {
                nvectors.push_back(std::vector<T>(component_size, default_));
            }
            npvectors = initial_npvectors;
            for (int i = 0; i < initial_npvectors; i++) {
                pvectors.push_back(std::vector<T>(component_size, default_));
            }
        }
        int nsize () {
            return nnvectors * component_size;
        }
        int psize () {
            return npvectors * component_size;
        }
        T& operator[] (int index) {
            if (index > max_pindex) {
                max_pindex = index;
            }
            else if (index < max_nindex) {
                max_nindex = index;
            }
            if (index >= 0) {
                if (index >= (*this).psize()) {
                    int x = (index - (*this).psize()) / component_size;
                    npvectors += x + 1;
                    for (int i = 0; i <= x; i++) {
                        pvectors.push_back(std::vector<T>(component_size, default_));
                    }
                }
                auto it = pvectors.begin();
                int j;
                for (j = 0; j < index / component_size; j++) {
                    ++it;
                }
                return (*it)[index - j * component_size];
            }
            else {
                index = -1 - index;
                if (index >= (*this).nsize()) {
                    int x = (index - (*this).nsize()) / component_size;
                    nnvectors += x + 1;
                    for (int i = 0; i <= x; i++) {
                        nvectors.push_back(std::vector<T>(component_size, default_));
                    }
                }
                auto it = nvectors.begin();
                int j;
                for (j = 0; j < index / component_size; j++) {
                    ++it;
                }
                return (*it)[index - j * component_size];
            }
        }
        int get_nnvectors() {
            return nnvectors;
        }
        int get_npvectors() {
            return npvectors;
        }
        int get_max_nindex() {
            return max_nindex;
        }
        int get_max_pindex() {
            return max_pindex;
        }
        int get_component_size() {
            return component_size;
        }
};