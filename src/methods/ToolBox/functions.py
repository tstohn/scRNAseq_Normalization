def removesuffix(word, suffix):
    if(word.endswith(suffix)):
        wordNew = word[:-len(suffix)]
        word = wordNew
    return(word)