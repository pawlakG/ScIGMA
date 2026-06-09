import os
import re
import json

french_words = [
    "pour", "avec", "dans", "données", "mise à jour", "ajout", "erreur",
    "récupérer", "sélection", "afficher", "ajouter", "supprimer", "mise", 
    "fichier", "utilisateur", "générer", "récupération", "déclencheur",
    "rafraîchit", "assignation", "immédiate", "verrouillée", "stricte",
    "désactivé", "déjà", "présent", "onglet", "actif", "bouton", "cliqué",
    "aucun", "nouveau", "nouvelle", "choix", "écrasement", "étape", "génération"
]
pattern = re.compile(r'\b(' + '|'.join(french_words) + r')\b', re.IGNORECASE)

results = {}

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            path = os.path.join(root, f)
            with open(path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
            
            file_results = []
            for i, line in enumerate(lines):
                if pattern.search(line) and not line.strip().startswith('library('):
                    file_results.append({"line": i+1, "content": line.strip('\n')})
            
            if file_results:
                results[path] = file_results

with open("french_lines.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
