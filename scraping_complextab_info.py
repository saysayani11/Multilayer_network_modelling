import requests
import os
from bs4 import BeautifulSoup

# Define the URL of the webpage containing the TSV files
url = 'http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/'

# Send a request to the webpage and get the HTML content
response = requests.get(url)
html = response.content

# Parse the HTML content with BeautifulSoup
soup = BeautifulSoup(html, 'html.parser')

# Find all the links on the webpage and extract their URLs
links = soup.find_all('a')
urls = [link.get('href') for link in links]

# Filter the URLs to only include TSV files
tsv_urls = [u for u in urls if u.endswith('.tsv')]

# Download each TSV file
for tsv_url in tsv_urls:
    # Construct the full URL of the TSV file
    full_url = url + tsv_url
    
    # Define the filepath to save the downloaded file
    filepath = os.path.join(r"C:\Users\saySa\OneDrive\Desktop\ftp-data", tsv_url)
    
    # Send a request to download the file and save it to disk
    response = requests.get(full_url)
    with open(filepath, 'wb') as f:
        f.write(response.content)
