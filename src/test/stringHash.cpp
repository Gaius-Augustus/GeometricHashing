#include <iostream>
#include <string>

int main(int argc, char * argv[]) {
    if (argc != 2) {return -1; }
    std::string s{argv[1]};
    std::hash<std::string> hasher;
    std::cout << hasher(s) << std::endl;
    return 0;
}
