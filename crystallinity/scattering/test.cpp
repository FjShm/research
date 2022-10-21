#include <iostream>
#include <vector>
#include <algorithm>

int main(void)
{
        std::vector<int> vec;
        vec.push_back(5);
        vec.push_back(2);
        vec.push_back(4);

        std::vector<int>::iterator itr;
        const int wanted = 1; // 2に対応するindex(この例では1)を取って来たい
        itr = std::find(vec.begin(), vec.end(), wanted);
        if (itr == vec.end()) std::cout << "search failed" << std::endl;
        const int wanted_index = std::distance(vec.begin(), itr);
        std::cout << "Result: vec[" << wanted_index << "] = " << *itr << std::endl;
}
