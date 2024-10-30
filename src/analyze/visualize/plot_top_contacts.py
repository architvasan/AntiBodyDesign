import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

resid_data = pd.read_csv(f'out/meta/contacts_meta0/prob_nmnat2_disord.dat')
# create DataFrame

print(resid_data[' contact prob'])
df = pd.DataFrame({'Metastate 0':list(resid_data[' contact prob'])},
                    index=resid_data['resid #'])



for it in range(1, 6):
    resid_data = pd.read_csv(f'out/meta/contacts_meta{it}/prob_nmnat2_disord.dat')
    df[f'Metastate {it}'] = list(resid_data[' contact prob'])


print(df)

# Plot the DataFrame
ax = df.plot(kind='bar', stacked=True, color=['red', 'skyblue', 'green', 'orange', 'magenta', 'brown'], figsize=(30, 10))

# Increase width of each bar
for container in ax.containers:
    for bar in container:
        bar.set_width(0.9)  # Adjust the width as needed


# Get the current x-tick positions
xticks = ax.get_xticks()

# Filter x-tick positions to show every 5th one
step = 5
filtered_xticks = xticks[::step]

# Set new x-ticks and format labels
ax.set_xticks(filtered_xticks)
ax.set_xticklabels([f'{int(resid_data["resid #"][int(tick)])}' for tick in filtered_xticks])#[f'{int(tick):,}' for tick in filtered_xticks])  # Format labels with commas
    
# Add title and labels
plt.xlabel('Resids', fontsize=14)
plt.ylabel('Contact Probability', fontsize=14)
plt.xticks(rotation=45)

# Show the plot
plt.savefig(f'images/contacts_nmnat2_disord_metas.png', bbox_inches='tight', dpi=300)
plt.close()



