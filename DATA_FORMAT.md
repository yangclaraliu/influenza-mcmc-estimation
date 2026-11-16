# Data format

The source workbooks are not included because they contain individual-level
health and demographic information. Place authorized local copies beside the
Python script using these exact filenames:

- `流感数据年龄分布.xlsx`
- `流感疫苗数据库3_小学.xlsx`
- `流感疫苗数据库4_初中.xlsx`
- `流感疫苗数据库5_高中.xlsx`

## School workbooks

Each workbook must contain a header row and at least eight columns. The script
uses these zero-based column positions:

- Position 1: influenza status; values greater than zero are treated as infected.
- Position 2: vaccination status; values greater than zero are treated as vaccinated.
- Position 5: school identifier.
- Position 7: grade number, starting at 1.

Records with missing or nonpositive school or grade values are removed. A
school is included only if it has records for every grade represented in that
workbook.

## Age-distribution workbook

The workbook must have at least 148 data rows. The script reads zero-based
position 8 for elementary schools, 9 for middle schools, and 10 for high
schools.

Values greater than 1 are interpreted as percentages and divided by 100.

Do not commit the original workbooks unless you have confirmed that sharing
them complies with participant consent, institutional policy, and applicable
privacy law.
