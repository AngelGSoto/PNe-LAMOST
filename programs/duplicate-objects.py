from astropy.coordinates import SkyCoord
from astropy.table import Table
import pandas as pd

# Sample data with RA and DEC values
dfTrue = pd.read_csv("DR7lamost-Likely-PN.csv")
table = Table.from_pandas(dfTrue)

# Convert RA and DEC values to SkyCoord objects
coordinates = SkyCoord(table['RA'], table['DEC'], unit='deg', frame='icrs')

# Create an empty list to store duplicate indices
duplicate_indices = []

# Compare each pair of objects to find duplicates
for i in range(len(coordinates)):
    for j in range(i + 1, len(coordinates)):
        if coordinates[i].separation(coordinates[j]).arcsecond < 1.0:  # Adjust the threshold as needed
            duplicate_indices.append((i, j))

# Print the duplicate object pairs
for i, j in duplicate_indices:
    print(f'Duplicate Objects: Object {table["FileName"][i]} and Object {table["FileName"][j]}')
    

# Create a new table without duplicates
new_dfTrue = table.copy()
new_dfTrue.remove_rows([index[1] for index in duplicate_indices])

# save
# Save the table to a CSV file
new_dfTrue.write('DR7lamost-Likely-PN-duplicate.ecsv', format='ascii.ecsv', overwrite=True)

