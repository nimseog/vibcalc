import requests

def get_url_text(url):
    headers = {'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:55.0) Gecko/20100101 Firefox/55.0'}
    return requests.get(url, headers=headers).text
