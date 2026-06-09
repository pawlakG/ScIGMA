import os

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        # String translations
        line = line.replace("Erreur : L'Assay 'compass_imputed' est introuvable. Veuillez relancer l'inférence COMPASS.", "Error : 'compass_imputed' Assay not found. Please rerun COMPASS inference.")
        line = line.replace("Erreur : Certains variants sélectionnés n'existent pas dans la matrice.", "Error : Some selected variants do not exist in the matrix.")
        line = line.replace("Les matrices ne contiennent aucun nom de colonne (cellules).", "Matrices do not contain any column names (cells).")
        line = line.replace("Fichiers introuvables pour le préfixe : %s. L'inférence a-t-elle convergé ?", "Files not found for prefix: %s. Did the inference converge?")
        line = line.replace("Aucun singulet trouvé. Matrice inexploitable.", "No singlets found. Matrix unusable.")
        
        # Comments deletion
        lower_line = line.lower()
        if '#' in line:
            if "aucun pré-requis n'est rempli" in lower_line \
               or "l'interface reste propre sans message" in lower_line \
               or "empêche l'erreur seq.default" in lower_line \
               or "aiguillage dynamique" in lower_line \
               or "agrégation spatiale" in lower_line \
               or "jointure sur le nouveau format" in lower_line:
                if line[:line.find('#')].strip() == '':
                    continue
                else:
                    line = line[:line.find('#')].rstrip() + '\n'
                    
        new_lines.append(line)
        
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            process_file(os.path.join(root, f))

print("Phase 2 purge completed.")
