#include <test_system/test_system.hpp>
#include <kernel/util/stringify.hpp>
#include <kernel/util/exception.hpp>

#include <cstdlib>
#include <iostream>

using namespace FEAST;
using namespace TestSystem;

int main(int /*argc*/, char** /*argv*/)
{
    int result(EXIT_SUCCESS);
    size_t list_size(0);

    /*for (TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
            i != i_end ; )
    {
            ++i;
    }*/

    list_size = TestList::instance()->size();
    unsigned long iterator_index(1);
    for (TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
            i != i_end ; )
    {
        CONTEXT("When running test case '" + (*i)->id() + "':");
        try
        {
            std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id() + " [Backend: "
                << (*i)->get_tag_name() << "]" << " [Precision: "<< (*i)->get_prec_name() << "]" << std::endl;
            (*i)->run();
            std::cout << "PASSED" << std::endl;
        }
        catch (TestFailedException & e)
        {
            std::cout << "FAILED: " << (*i)->id() << std::endl << stringify(e.what()) << std::endl;
            result = EXIT_FAILURE;
        }
        i = TestList::instance()->erase(i);
        iterator_index++;
    }

    return result;
}
