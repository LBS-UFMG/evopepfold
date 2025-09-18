import requests

def runAlphaFold(fasta_path):
    # The URL of the Flask API (replace <YOUR_SERVER_IP> with the actual server IP)
    url = 'http://localhost:5000/run_alphafold'

    # Path to the fasta file
    fasta_file_path = 'inputs.fasta'

    # Send the POST request with the fasta file
    with open(fasta_file_path, 'rb') as fasta_file:
        files = {'fasta_file': fasta_file}
        response = requests.post(url, files=files)

    # Print the response
    print(response.json())

def runAlphaFold_fake(fasta_path):
    pass