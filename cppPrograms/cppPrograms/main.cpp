#include <iostream>
#include <string>
#include "hello.hpp"

using namespace std;

void sample(string m)
{
    cout << m + " rain";
    return;
}
int main()
{
    string entry;
    getline(cin, entry);
    helloWorld();
    sample(entry);
    return 0;
}
