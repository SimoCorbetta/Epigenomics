#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from adjustText import adjust_text


# In[2]:


input_file = sys.argv[1]
output_file = sys.argv[2]

# In[3]:


with open(input_file, 'r') as f:
    x = f.readlines()
    lista_di_liste = [linea.strip().split('\t') for linea in x]
    states = []
    for l in lista_di_liste:
        #state = l[9]
        state = l[13]
        states.append(state)


# In[4]:


states


# In[5]:


term_count = {term: states.count(term) for term in set(states)}
term_count


# In[6]:


# Ordinamento dei valori in ordine crescente
sorted_dict = dict(sorted(term_count.items(), key=lambda item: item[1]))
sorted_dict


# In[7]:


# Dati forniti
data_new = [
    ("TssA", "Active TSS", "Red", "255,0,0"),
    ("TssFlnk", "Flanking TSS", "Orange Red", "255,69,0"),
    ("TssFlnkU", "Flanking TSS Upstream", "Orange Red", "255,69,0"),
    ("TssFlnkD", "Flanking TSS Downstream", "Orange Red", "255,69,0"),
    ("Tx", "Strong transcription", "Green", "0,128,0"),
    ("TxWk", "Weak transcription", "DarkGreen", "0,100,0"),
    ("EnhG1", "Genic enhancer1", "GreenYellow", "194,225,5"),
    ("EnhG2", "Genic enhancer2", "GreenYellow", "194,225,5"),
    ("EnhA1", "Active Enhancer 1", "Orange", "255,195,77"),
    ("EnhA2", "Active Enhancer 2", "Orange", "255,195,77"),
    ("EnhWk", "Weak Enhancer", "Yellow", "255,255,0"),
    ("ZNF/Rpts", "ZNF genes & repeats", "Medium Aquamarine", "102,205,170"),
    ("Het", "Heterochromatin", "PaleTurquoise", "138,145,208"),
    ("TssBiv", "Bivalent/Poised TSS", "IndianRed", "205,92,92"),
    ("EnhBiv", "Bivalent Enhancer", "DarkKhaki", "189,183,107"),
    ("ReprPC", "Repressed PolyComb", "Silver", "128,128,128"),
    ("ReprPCWk", "Weak Repressed PolyComb", "Gainsboro", "192,192,192"),
    ("Quies", "Quiescent/Low", "White", "240,240,240")
]

# Creazione del dizionario con la prima colonna come chiave e valori RGB normalizzati
color_dict_new_normalized = {row[0]: tuple(int(val) / 255 for val in row[3].split(',')) for row in data_new}


# In[8]:


terms = list(term_count.keys())
counts = list(term_count.values())

colors = [color_dict_new_normalized[term] for term in terms]
plt.figure(figsize=(8, 7))
plt.bar(terms, counts, color=colors)
plt.xlabel('Chrm states')
plt.ylabel('Count')
plt.title('States in my peaks')
plt.xticks(terms, rotation=45, ha='right')
plt.savefig(f'states_barplot_{output_file}.pdf', format='pdf') # to save the pdf of the peaks
#plt.show()


# In[12]:


terms = list(sorted_dict.keys())
counts = list(sorted_dict.values())

min_percentage = 1.5

# Generazione del grafico a torta
fig, ax = plt.subplots(figsize=(8.5, 6))
wedges, _ = ax.pie(counts, colors=[color_dict_new_normalized[term] for term in terms], startangle=90,
                   wedgeprops=dict(width=0.5, edgecolor='w'))

# Calcolo delle percentuali
total = sum(counts)
percentages = [100 * count / total for count in counts]

# Calcolo delle percentuali
for i, (wedge, percentage) in enumerate(zip(wedges, percentages)):
    if percentage > min_percentage:
        ang = (wedge.theta2 - wedge.theta1) / 2.0 + wedge.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontal_alignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connection_style = f"angle,angleA=0,angleB={ang}"
        ax.annotate(f"{terms[i]}: {percentage:.1f}%", (x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontal_alignment, fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black"),
                    arrowprops=dict(arrowstyle="-", connectionstyle=connection_style))

# Stati con percentuale inferiore al 1%
small_percentage_data = [(term, percentage) for term, percentage in zip(terms, percentages) if percentage <= min_percentage]

# Ordinamento degli stati in base alla percentuale in ordine decrescente
small_percentage_data_sorted = sorted(small_percentage_data, key=lambda x: x[1], reverse=True)

# Creazione degli elementi della legenda per gli stati con percentuale minore di 1% in ordine decrescente
legend_patches = [mpatches.Patch(color=color_dict_new_normalized[term], label=f'{term}: {percentage:.1f}%')
                  for term, percentage in small_percentage_data_sorted]
# Plot della legenda per gli stati con percentuale minore di 1%
if len(legend_patches) > 0:
    ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1, 1)) 

#plt.title('Distribution of States in Peaks')
plt.savefig(f'states_pie_chart_{output_file}.pdf', format='pdf')  # Salva l'immagine del grafico a torta
#plt.show()

