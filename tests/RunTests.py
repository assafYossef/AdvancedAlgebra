from alive_progress import alive_bar

from PrimeFieldTests import PrimeFieldTests
from FiniteFieldTests import FiniteFieldTests

# This file will contain all the tests for the classes


def generic_run(public_method_names, p):
    failed_tests = []
    with alive_bar(len(public_method_names), spinner='squares') as bar:  # declare your expected total
        # run all tests from class
        for method in public_method_names:
            try:
                print(f"Running {method}")
                getattr(p, method)()
                bar()
            except Exception as e:
                failed_tests.append(method)
                print(str(e) + " on Test: " + method)

    print(f"{len(public_method_names) - len(failed_tests)} Succeed, {len(failed_tests)} Failed")
    if (len(failed_tests) > 0):
        print("The tests that failed are:")
        print(failed_tests)



def run_all_tests():
    p = PrimeFieldTests(p = 11)
    public_method_names = [method for method in dir(p) if callable(getattr(p, method)) if method.startswith('test')]
    generic_run(public_method_names, p)

    p = FiniteFieldTests()
    public_method_names = [method for method in dir(p) if callable(getattr(p, method)) if method.startswith('test')]
    generic_run(public_method_names, p)


if __name__ == "__main__":
    run_all_tests()

