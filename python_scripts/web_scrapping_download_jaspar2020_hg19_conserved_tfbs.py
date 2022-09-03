from bs4 import BeautifulSoup
import requests
import os
page = requests.get("http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/hg19/")
#print(page.content)

soup = BeautifulSoup(page.content, 'html.parser')
#print(soup.prettify())#print formated html data
#print(soup.get_text())#print text from the html page
#print(soup.find_all('p'))
#print(soup.find_all('p', class_='outer-text'))
#print(soup.find_all(class_="outer-text"))
#print(soup.find_all(id="first"))

#Write html file in text:
with open("output1.txt", "w") as file:
    file.write(str(soup))

for a in soup.find_all('a', href=True):
    if(a['href'][0] == 'M'):
        cmd = 'wget expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/hg19/'+a['href']
        os.system(cmd)
