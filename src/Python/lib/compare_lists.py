import sys
import os


def main():

    if len(sys.argv) == 1:
        print(f"Missing required required arguments. \nusage: py compare_lists.py <file one> <file two>")
        return

    if len(sys.argv) == 2:
        print(f"Missing a second argument for the second file name.\nusage: py compare_lists.py <file one> <file two>")
        return
    
    lines1 = set()
    lines2 = set()

    if not os.path.exists(sys.argv[1]):
        print(f"File Argument 1 <<{sys.argv[1]}>> does not exist.")
        return

    if not os.path.exists(sys.argv[2]):
        print(f"File Argument 2 <<{sys.argv[2]}>> does not exist.")
        return

    with open(sys.argv[1]) as f1:
        with open(sys.argv[2]) as f2:
            lines1 = {x.strip() for x in set(f1.readlines())}
            lines2 = {x.strip() for x in set(f2.readlines())}

    print(f"# lines in file 1: {len(lines1)}")
    print(f"# lines in file 2: {len(lines2)}")

    repeated = {x for x in lines1 if x in lines2}

    print(f"Num repeated: {len(repeated)}")
    # print(f"Repeated: {repeated}")

    print(f"Lines in file 1, but not in file 2:\n{lines1-lines2}")
    print(f"Lines in file 2, but not in file 1:\n{lines2-lines1}")

if __name__ == "__main__":
    main()
