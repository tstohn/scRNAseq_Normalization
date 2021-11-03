import sys

def removesuffix(word, suffix):
    if(word.endswith(suffix)):
        wordNew = word[:-len(suffix)]
        word = wordNew
    return(word)

def printToTerminalOnce(string):
    tmp = sys.stdout
    sys.stdout = sys.__stdout__
    print(string)
    sys.stdout = tmp