import os
import re

french_words = [
    "pour", "avec", "dans", "données", "mise à jour", "ajout", "erreur",
    "récupérer", "sélection", "afficher", "ajouter", "supprimer", "mise", 
    "fichier", "utilisateur", "générer", "récupération", "déclencheur",
    "rafraîchit", "assignation", "immédiate", "verrouillée", "stricte"
]
pattern = re.compile(r'\b(' + '|'.join(french_words) + r')\b', re.IGNORECASE)

results = {"comments": 0, "strings": 0, "files": set()}

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            path = os.path.join(root, f)
            with open(path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
                for line in lines:
                    if pattern.search(line):
                        results["files"].add(path)
                        # Check if comment
                        if '#' in line:
                            # Simplistic check: is the match after #?
                            comment_part = line[line.find('#'):]
                            string_part = line[:line.find('#')]
                            if pattern.search(comment_part):
                                results["comments"] += 1
                            if pattern.search(string_part) and ('"' in string_part or "'" in string_part):
                                results["strings"] += 1
                        elif '"' in line or "'" in line:
                            results["strings"] += 1

print(f"Total files with French: {len(results['files'])}")
print(f"Total French comments: {results['comments']}")
print(f"Total French UI/Strings: {results['strings']}")
