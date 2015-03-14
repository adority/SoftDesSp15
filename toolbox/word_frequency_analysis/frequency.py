""" Analyzes the word frequencies in a book downloaded from
	Project Gutenberg """

import string

def get_word_list(file_name):
	""" Reads the specified project Gutenberg book.  Header comments,
		punctuation, and whitespace are stripped away.  The function
		returns a list of the words used in the book as a list.
		All words are converted to lower case.
	"""
	f = open(file_name,'r')
	lines = f.readlines()
	curr_line = 0
	while lines[curr_line].find('CHAPTER 1') == -1:
		curr_line += 1
	lines = lines[curr_line+1:]
	f.close()

	words = []
	for line in lines:
		for word in line.split():
			word = word.strip(string.punctuation)
			word = word.lower()
			if word != '':
				words.append(word)

	return words

#get_word_list('20000leagues.txt')


def get_top_n_words(word_list,n):
	""" Takes a list of words as input and returns a list of the n most frequently
		occurring words ordered from most to least frequently occurring.

		word_list: a list of words (assumed to all be in lower case with no
					punctuation
		n: the number of words to return
		returns: a list of n most frequently occurring words ordered from most
				 frequently to least frequentlyoccurring
	"""
	frequencies = {}
	for word in word_list:
		if word in frequencies:
			frequencies[word] += 1
		else:
			frequencies[word] = 1
	ordered_by_frequency = sorted(frequencies, key=frequencies.get, reverse=True)
	return ordered_by_frequency[0:n]

word_list = get_word_list('20000leagues.txt')
print get_top_n_words(word_list,100)

