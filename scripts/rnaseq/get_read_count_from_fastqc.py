from bs4 import BeautifulSoup
import sys

def main(fn):
    with open(fn, 'rb') as f:
        soup = BeautifulSoup(f.read(), 'html.parser')
    t = soup.find('td', text='Total Sequences')
    arr = t.fetchNextSiblings()
    if len(arr) != 1:
        raise AttributeError("Expected to find a single sibling")
    return int(arr[0].text)


if __name__ == "__main__":
    fn = sys.argv[1]
    n = main(fn)
    print "%s, %d" % (fn, n)