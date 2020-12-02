import camelot
tables = camelot.read_pdf('/Users/biplab/Desktop/H-Canada-PBiT.pdf')
tables.export('foo.csv', f='csv', compress=True)
print(tables[0])