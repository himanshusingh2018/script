def extract_pmid_by_keyword(search_term, start_date, nsearch):
    '''
    Input:
        search_term = "Alzheimer"
        nsearch = 150; number of pubmed id to be searched
    Output:
        list of PMID
    RUN:
        extract_pmid_by_keyword(search_term)
        Return list of PMIDs which has Alzheimer in the Abstract
    '''
    import re
    from Bio import Entrez
    Entrez.email = "himanshu720@gmail.com"# Entrez email (provide your email address)
    handle = Entrez.esearch(db="pubmed", term=search_term, mindate=start_date, retmax=nsearch)# Search PubMed with specified term and retrieve the PubMed IDs
    record = Entrez.read(handle)
    pmids = record["IdList"]
    return(pmids)

def extract_email_from_url(url):
    import re
    import requests

    try:
        response = requests.get(url)  # Send a GET request to the web page
        response.raise_for_status()  # Raise an exception if the response status is not successful
        email_pattern = r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b"  # Extract email addresses using regular expressions
        emails = list(set(re.findall(email_pattern, response.text)))
        return emails
    except requests.exceptions.RequestException as e:
        print(f"Error occurred while accessing the URL: {url}")
        print(e)
        return []

def extract_url_from_pmid(pmid):
    '''
    Input:
        pmid: PubMed ID
    Output:
        DOI (Digital Object Identifier)
    '''
    import requests
    from Bio import Entrez
    
    Entrez.email = "himanshu720@gmail.com"  # Provide your email address for Entrez
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")  # Fetch the PubMed record in XML format
    record = Entrez.read(handle)["PubmedArticle"][0]
    article = record["PubmedData"]["ArticleIdList"]
    doi = None
    # Iterate through the article identifiers to find the DOI
    for identifier in article:
        if identifier.attributes["IdType"] == "doi":
            doi = identifier
    return(f'https://doi.org/{doi}')

pmids = extract_pmid_by_keyword(search_term='Alzheimer', start_date=2015, nsearch=500)
print(len(pmids))

with open("output.txt", "w") as file:
    for pmid in pmids:
        url = extract_url_from_pmid(pmid)
        emails = extract_email_from_url(url)
        output = f"{pmid} {emails}\n"
        file.write(output)


