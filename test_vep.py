import urllib.request
import json
url = "https://rest.ensembl.org/vep/human/hgvs?phenotypes=1"
data = json.dumps({"hgvs_notations": ["7:g.140453136A>T"]}).encode('utf-8')
req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"})
response = urllib.request.urlopen(req)
result = json.loads(response.read().decode('utf-8'))
print(json.dumps(result[0].get('colocated_variants', []), indent=2))
