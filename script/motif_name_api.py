import pandas as pd
import requests
import os


def convert_motif_id_to_tf_name(motif_id):
    # JASPAR REST API endpoint for motif information
    api_url = 'https://jaspar.genereg.net/api/v1/matrix/{}/'.format(motif_id)

    try:
        # Make a GET request to the API
        response = requests.get(api_url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the JSON response
            data = response.json()

            # Extract the transcription factor name
            tf_name = data['name']

            return tf_name
        else:
            # Handle error cases
            print("Error: {}".format(response.status_code))
            return None
    except Exception as e:
        print("An error occurred: {}".format(e))
        return None

# List of input files
input_files = ['GPC_RG_chromvar.tsv', 'astro_opc_chromvar.tsv', 'astro_GPC_chromvar.tsv', 'opc_GPC_chromvar.tsv']

for input_file in input_files:
    # Read the input table with row names
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    # Add a new column for transcription factor names using the row names (index)
    df['TF_name'] = df.index.map(convert_motif_id_to_tf_name)

    # Save the updated table
    output_file = os.path.splitext(input_file)[0] + '_with_TF_names.tsv'
    df.to_csv(output_file, sep='\t')

