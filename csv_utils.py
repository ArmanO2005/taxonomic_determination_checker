import pandas as pd
import re
import string
from kew_loader_utils import *

def update_determinations(path_to_csv, tree, genus_error_dist=3, species_error_dist=3):
    """
    this function takes the file path to a CSV with one column 

    The column should be called "name" and include all names to be checked/updated

    the function returns a csv with three columns, 
    the original names,
    the spell checked names,
    the accepted names from kew if the spell checked name is a synonym,
    and the author
    """
    df = pd.read_csv(path_to_csv)

    df = TaxonNoAuthor(df, 'name')

    df['checked_name'] = df['***no_author'].apply(lambda x: try_helper(tree.query, x.strip().lower(), genus_error_dist, species_error_dist))

    df[['accepted_name', 'author']] = df['checked_name'].apply(lambda x: pd.Series(tree.getAcceptedName(x) if x else (None, None)))

    df = df.drop(columns=['binomial_match', 'additions_match', '***no_author'])

    df.to_csv('output.csv', index=False)


def TaxonNoAuthor(data, columnName):
    """
    Modifies the dataframe data to contain one new columns:
    ***no_author: the binomial with the author name ommitted
    """
    binomialPattern = (
        r'([A-Z][a-z]+ \'.+\'|[A-Z][a-z]+ aff\. [a-z-]+|'
        r'[A-Z][a-z]+ cf\. [a-z-]+|[A-Z][a-z]+ [a-z-]+)'
    )
    
    additionsPattern = (
        r'(ssp\. [a-z]+|subsp\. [a-z]+|var\. [a-z]+|\'[a-z]+\'|\'.+\')'
    )

    def get_first_match(pattern, text):
        matches = re.findall(pattern, text)
        return matches[0] if matches else None

    def get_all_matches(pattern, text):
        matches = re.findall(pattern, text)
        return matches if matches else []
    
    def strip_sp(text):
        pattern = r"(?:\s[Ss][Pp]\s|\s[Ss][Pp]\.)"
        if ' sp' in text or ' Sp' in text:
            edited = re.sub(pattern, ' ', text)
            edited = edited.replace('  ', ' ').strip()
            return edited
        return text

    data["binomial_match"] = data[columnName].apply(lambda x: get_first_match(binomialPattern, str(x).strip().replace('  ', ' ')))
    data["additions_match"] = data[columnName].apply(lambda x: get_all_matches(additionsPattern, str(x).strip().replace('  ', ' ')))
    data["binomial_match"] = data["binomial_match"].fillna(data[columnName])
    
    for i in range(len(data)):
        if len(data.at[i, "additions_match"]) > 0 and len(data.at[i, "binomial_match"].split(" ")) > 1:
            if data.at[i, "additions_match"][0] == data.at[i, "binomial_match"].split(" ", 1)[1]:
                data.at[i, "additions_match"] = data.at[i, "additions_match"][1:]

    data["***no_author"] = data["binomial_match"] + data["additions_match"].apply(lambda x : ' ' + ' '.join(x))
    data["***no_author"] = data["***no_author"].apply(lambda x: str(x).translate(string.punctuation))
    data["***no_author"] = data["***no_author"].apply(lambda x: strip_sp(x).replace('  ', ' ').strip())

    return data

def try_helper(fxn, *args):
    try:
        return fxn(*args)
    except:
        return None