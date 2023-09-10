import requests
import os
os.chdir(r"C:\Users\saySa\OneDrive\Desktop\task_1\MSA\consurf_data")

# List of URLs you provided
urls = [         
    "https://consurfdb.tau.ac.il/DB/6K8KG/6K8KG_consurf_summary.txt"
 ]

# Directory where you want to save the files
save_directory = "local_folder"

# Create the directory if it doesn't exist
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# Loop through each URL
for url in urls:
    response = requests.get(url)
    if response.status_code == 200:
        data = response.text
        
        # Create a filename based on the URL. Here, I'm using the last part of the URL.
        filename = os.path.join(save_directory, url.split("/")[-1])
        
        # Save the data to a local file
        with open(filename, 'w') as file:
            file.write(data)
        
        print(f"Data from {url} saved to {filename}")
    else:
        print(f"Failed to retrieve data from {url}. Status code: {response.status_code}")
