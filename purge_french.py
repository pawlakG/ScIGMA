import os
import re

french_words = [
    " pour ", " avec ", " dans ", "données", "mise à jour", " erreur",
    "récupérer", "sélection", "afficher", "ajouter", "supprimer", "fichier",
    "utilisateur", "générer", "récupération", "déclencheur", "rafraîchit",
    "assignation", "immédiate", "verrouillée", "stricte", "désactivé",
    " déjà", " présent", " onglet", " actif ", "bouton", "cliqué", "écrasement",
    "étape", "génération", "sécurité", "éviter", "mise à l'échelle",
    "cohérence", "fermeture", "balises", "saturer", "sauvegardé", "préparation",
    "formatage", "enfants", "traités", "dernier", "écraser", "garantir",
    "couleurs", "discrètes", "pivotage", "transposition", "obligatoire",
    "itérer", "lignes", "retourne", "matrice", "restaurer", "création", "isolée",
    "correspondre", "celui", "défini", "voir ", "événement", "déclenché",
    "recacher", " avis ", "recharge", "simples", "absente", "assurer",
    "niveaux", "transfère", "partagé", "préserve", "référence", "mémoire",
    "injection", "nouvelles", "vérification", "annoté", "atomisation",
    "éviter", "anciens", "brutes", "bruts", "révèle", " force ", "désactivé",
    "récupérer", "écoute", "nouvel", "interrupteur", "rafraîchit", "termine",
    " actif", "pare-feu", "imputées", "invalidation", "re-générés", "prendre",
    "largeur", "garantit", "création", "compatible", "exécution", "lancement",
    "réel ", " valide", "stoppe", " cliqué", "début", "existant", " fin ",
    "expulse", "doublets", "purge", "autoriser", "stricte",
    "fusion", "petits", "conservation", "mise en "
]

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        # String translations
        line = line.replace('Aucune sous-population créée à enregistrer.', 'No sub-population created to save.')
        line = line.replace('Sélection : ', 'Selection: ')
        line = line.replace('1/2 - Extraction HDF5 et préparation des matrices...', '1/2 - Extracting HDF5 and preparing matrices...')
        line = line.replace('Échec critique C++ :', 'Critical C++ failure:')
        
        # Detect French comment
        is_french_comment = False
        if '#' in line:
            comment_part = line[line.find('#'):].lower()
            for w in french_words:
                if w in comment_part:
                    is_french_comment = True
                    break
                    
        if is_french_comment:
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

print("Purge completed.")
