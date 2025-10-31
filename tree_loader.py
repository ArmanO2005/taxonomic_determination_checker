import pandas as pd
from BK_Tree import BKTree
import pickle
import re
import string
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class kew_tree:
    def __init__(self, kew_data_df):
        self.kew_data = kew_data_df
        self.genus_tree = BKTree()
        self.mapper_dict = {}        

    def build(self, prune=None):
        """
        prune: type(pandas series of binomial names), if specified, will only build the tree for genera in the prune list, 
        useful if you're only checking one list of names and dont plan on reusing the data structure
        """
        self.kew_data['genus'] = self.kew_data['genus'].str.strip().str.lower()
        self.kew_data['specificepithet'] = self.kew_data['specificepithet'].str.strip().str.lower()

        genera_list = self.kew_data['genus'].unique().tolist()

        if prune is not None:
            prune = pd.DataFrame({'name': prune}).reset_index(drop=True)
            prune = TaxonNoAuthor(prune, 'name')
            prune['genus'] = prune["***no_author"].apply(lambda x: x.split()[0].lower())
            prune_genera = prune['genus'].unique().tolist()
            # genera_list = [genus for genus in genera_list if genus in prune_genera]
            # print(len(genera_list) / len(prune_genera))

        
        for genus in genera_list:
            self.genus_tree.insert(str(genus))

        if prune is not None:
            checked_prune_genera = []
            for genus in prune_genera:
                matches = self.genus_tree.search(str(genus), 2, early_termination=False)
                if matches:
                    matched_genus = [m[0] for m in matches]
                    checked_prune_genera += matched_genus
            genera_list = pd.Series(checked_prune_genera).unique().tolist()

        for genus in genera_list:
            genus_df = self.kew_data[self.kew_data['genus'] == genus]
            species_list = genus_df['scientfiicname'].apply(lambda x: " ".join(x.split()[1:]))
            species_list = species_list.unique().tolist()

            species_tree = BKTree()
            for species in species_list:
                species_tree.insert(str(species))

            self.mapper_dict[genus] = species_tree

    def query(self, formatted_binomial, error_dist_genus, error_dist_species, early_termination=False):
        """
        binomial must be stripped of author stuff. Must start with Genus and then species and then subspecific terms. If its a cultivar w/o a specific name i.e. Rhododendron "CULTIVAR NAME" or has a non-standrd
        naming convention, it wont work
        """
        genus = formatted_binomial.split()[0]

        matched_genus = self.genus_tree.search(genus, error_dist_genus, early_termination=early_termination)
        if not matched_genus:
            logger.warning(f"No matching genus found for '{genus}'. Adjust the error distance or check the spelling.")
            return None

        matched_genus = sortOutput(matched_genus)
        matched_genus = matched_genus[0][0]

        try:
            species = " ".join(formatted_binomial.split()[1:])
        except:
            return matched_genus

        species_tree = self.mapper_dict[matched_genus]
        matched_species = species_tree.search(species, error_dist_species, early_termination=early_termination)
        if not matched_species:
            return matched_genus
        

        matched_species = sortOutput(matched_species)
        matched_species = matched_species[0][0]
        return matched_genus + " " + matched_species
    
    def getAcceptedName(self, name):
        name = name.capitalize().strip()
        row = self.kew_data[self.kew_data['scientfiicname'] == name]
        if len(row) > 1:
            row = row[row['taxonomicstatus'] == 'accepted']
        if row.empty:
            logger.warning(f"Accepted name for '{name}' not found")
            return None, None, None, None
        
        if not row.empty:
            accepted_id = row.iloc[0]['acceptednameusageid']
            checked_synonym = row.iloc[0]["scientfiicname"]

            accepted_row = self.kew_data[self.kew_data['taxonid'] == accepted_id].iloc[0]

            if accepted_row['family']:
                family = accepted_row['family']
            else:
                family = None

            if accepted_row['scientfiicnameauthorship']:
                author = accepted_row['scientfiicnameauthorship']
            else:
                author = None

            if accepted_row['scientfiicname']:
                accepted_name = accepted_row['scientfiicname']
            else:
                accepted_name = None

            return checked_synonym, accepted_name, author, family
        


    # def save_tree(self):
    #     with open("built_object/kew_tree.pkl", 'wb') as f:
    #         pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    # def open_tree(self):
    #     with open("built_object/kew_tree.pkl", 'rb') as f:
    #         return pickle.load(f)


def sortOutput(list):
    if list:
        list.sort(key=lambda x: x[1])   
    return list

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


