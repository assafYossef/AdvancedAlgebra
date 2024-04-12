import os
import webbrowser
import subprocess
import time

from tests.RunTests import run_all_tests



def open_docs():
    docs_location = os.path.dirname(__file__) + "/docs"

    command = docs_location + "/make.bat"
    subprocess.check_call([command, "html"])

    index_html_file_location = docs_location + "/_build/html/index.html"

    webbrowser.open_new(index_html_file_location)


def about_us():
    print("Hello, World!")




def main():
    print("Welcome to advenced algebra final project! 😁")
    print("Created by Ido Nahum and Assaf Yossef")
    print("This project is about extended fileds and galois theory")
    print()
    print("Here you have a little user interface that we created for you")
    done = False
    while not done:
        print()
        print("Type '1' for PrimeFieldElement")
        print("Type '2' for FiniteField and FiniteFiledElement")
        print("Type '3' to run tests that we created 🧪")
        print("Type '4' to see the docs ℹ")
        print("Type '5' to see about us 👥")
        print("Type '6' to end")
        user_choise = input("Please enter your choise:")

        match user_choise:
            case "1":
                print(1)
            case "2":
                print(2)  
            case "3":
                run_all_tests()
            case "4":
                open_docs()
            case "5":
                about_us()
            case "6":
                done = True
            case _:
                print("Invalid input 😵, try again...")
        
        time.sleep(1)
    
    print("Bye Bye 👋! see you later")


if __name__ == "__main__":
    main()
