# from urllib.request import urlopen

# link = "http://www.vibrationdata.com/elcentro_NS.dat"

# f = urlopen(link)
# myfile = f.read()
# print(myfile)


# import requests

# link = "http://www.vibrationdata.com/elcentro_NS.dat"
# f = requests.get(link)
# print(f.text)

# import urllib2  # the lib that handles the url stuff

# data = urllib2.urlopen(target_url) # it's a file like object and works just like a file
# for line in data: # files are iterable
#     print(line)

target_url = 'http://www.vibrationdata.com/elcentro_NS.dat'
import urllib.request
data = urllib.request.urlopen(target_url)