from pattern.web import *
from facepy import exceptions
from facepy import GraphAPI
from pattern.en import sentiment, modality
f = Facebook(license = 'CAAEuAis8fUgBAH0onvYNvLzG7mBIRHqkQ6pqBFdqaafJPzzsXWsyJ9Eohm3UZBymwCwZBFDm3Q24kHEEVeQXhYeV0UgZCZAQLIerOBROZCyztZCViMwaw7FZC8CHHFRwgYumZBqGIx1KTw9nxZCQ8yiGrgWE9nZBqgcOqHinm0kwZCBeK54zpvW4NWMIe7pRZAHnZBgsZD')
token = 'CAACEdEose0cBAJdaCbPOAcpAF72QcVARtBTgaegxcESkVC7m3iWVxu0ahk9m5Fk6CA4912NRlrgojZAclIMQe1KsfZAS7qQWTzmPEwOFZBqqMQZCoJDbYcFCgqjIrdW1XIDiIQUoF7De7woDDIvLDuKBSiGeF4j0rTTkuGoaGEq6vOyz5s1gOue7Tm6YvGWqLZAkmK8VfzvW01saDPz3ZAgfsx3vBFBsoZD'
#This token will only last an hour and needs to be replaced with an access token found at https://developers.facebook.com/tools/explorer/
#A functional token can be obtained by clicking 'Get Access Token' and checking the box marked 'read_mailbox' under Extended Permissions
def messages():
	#sorts friends by sentiment and modality of their last message to you. Returns rankings as "Friends' Happiness" and "Friends' Confidence"
	graph = GraphAPI(token)
	me = f.profile()
	happiness = {}
	confidence = {}
	snippets = graph.fql('SELECT snippet, snippet_author FROM thread WHERE folder_id = 0 OR folder_id = 1 Limit 10000',3)
	#the above code was heavily influenced by arofcoding.blogspot.com/2012/10/python-script-to-fetch-messages-from.html
	#returns a dictionary of message snippets along with the corresponding facebook friend IDs
	for dictionary in snippets['data']:
	#puts snippets in a dictionary where each author is mapped to the sentiment of their message
		happiness[sentiment(dictionary['snippet'])] = dictionary['snippet_author']
		confidence[modality(dictionary['snippet'])] = dictionary['snippet_author']
	#ranks dictionary entries by positivity of sentiment
	happiness_rankings = rank(happiness)
	confidence_rankings = rank(confidence)
	print "Friends' Happiness (low to high):" 
	print happiness_rankings
	print "Friends' Confidence (low to high):"
	print confidence_rankings

def rank(dictionary):
#puts dictionary entries in order by lowest to highest key value
	rankings = []
	positivity = sorted(dictionary)
	for entry in positivity:
		rankings.append(dictionary[entry])
	return rankings

messages()