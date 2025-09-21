import pandas as pd
from BK_Tree import BKTree


class kew_tree:
    def __init__(self, kew_data_df):
        self.kew_data = kew_data_df
        self.genus_tree = BKTree()
        self.mapper_dict = {}        

    def build(self):
        self.kew_data['genus'] = self.kew_data['genus'].str.strip().str.lower()
        self.kew_data['specificepithet'] = self.kew_data['specificepithet'].str.strip().str.lower()

        genera_list = self.kew_data['genus'].unique().tolist()

        for genus in genera_list:
            self.genus_tree.insert(str(genus))

        for genus in genera_list:
            genus_df = self.kew_data[self.kew_data['genus'] == genus]
            species_list = genus_df['specificepithet'].unique().tolist()

            species_tree = BKTree()
            for species in species_list:
                species_tree.insert(str(species))

            self.mapper_dict[genus] = species_tree

    def query(self, formatted_binomial, error_dist_genus, error_dist_species):
        """
        binomial must be stripped of author stuff. Must start with Genus and then species and then subspecific terms. If its a cultivar or has a dumb, non-standrd
        naming convention, it wont work, im not a magician
        """
        genus = formatted_binomial.split()[0]

        matched_genus = self.genus_tree.search(genus, error_dist_genus)
        if not matched_genus:
            raise FileNotFoundError("No matching genus found. Adjust the error distance or check the spelling.")

        matched_genus = sortOutput(matched_genus)
        matched_genus = matched_genus[0][0]  # Take the closest match

        try:
            species = " ".join(formatted_binomial.split()[1:])
        except:
            return matched_genus, None

        species_tree = self.mapper_dict[matched_genus]
        matched_species = species_tree.search(species, error_dist_species)
        if not matched_species:
            return matched_genus, None
        

        matched_species = sortOutput(matched_species)
        matched_species = matched_species[0][0]  # Take the closest match
        return matched_genus, matched_species
    
    def getAcceptedName(self, name):
        if type(name) != str:
            name = " ".join(name).capitalize()
        row = self.kew_data[self.kew_data['scientfiicname'] == name]
        if not row.empty:
            accepted_id = row.iloc[0]['acceptednameusageid']

            accepted_row = self.kew_data[self.kew_data['taxonid'] == accepted_id].iloc[0]
            return accepted_row['scientfiicname'], accepted_row['scientfiicnameauthorship']

        return None


def sortOutput(list):
    if list:
        list.sort(key=lambda x: x[1])   
    return list