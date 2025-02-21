from os.path import join
import pandas as pd

DATA_MODELS = [
    "dataset",
    "sharingPlans",
    "education",
    "file",
    "grant",
    "person",
    "publication",
    "study",
    "tool",
]

COLS_TO_RENDER = [
    'Attribute',
    'Description',
    'Required',
    'Validation Rules',
    'Examples'
]

def on_pre_build(config, **kwargs) -> None:
    """Pre-process template files.
    
    Desired markdown: render model template so that
        - it is known which attributes require valid values
        - clicking on attribute will direct to valid values table
    """
    for model in DATA_MODELS:
        parent = join("modules", model)
        # Read both annotation properties and examples
        annotation_df = pd.read_csv(join(parent, 'annotationProperty.csv'), quoting=1).fillna("")
        examples_df = pd.read_csv(join(parent, 'exampleColumn.csv'), quoting=1).fillna("")

        # First select only the columns we want from annotation_df
        df = annotation_df[['Attribute', 'Description', 'Required', 'Validation Rules', 'Valid Values']]

        # Then add the Example column and rename it to Examples
        df = df.merge(
            examples_df[['Attribute', 'Example']], 
            on='Attribute', 
            how='left'
        ).rename(columns={'Example': 'Examples'})
        
        # If attribute has a list of valid values, create a link.
        for _, row in df[df['Valid Values'].ne("")].iterrows():
            attr_link = "[" + row['Attribute'] + (
                f"](../valid_values/{model}.md#attribute-"
                f"{row['Attribute'].lower().replace(' ', '-')})")
            df.at[_, 'Attribute'] = attr_link
        
        # For any validation rules with a regex, replace `\` with `\\`
        # for proper rendering.
        df['Validation Rules'] = (
            df['Validation Rules']
            .replace(r"\\", r"\\\\", regex=True))

        # Indicate "None" if there are no validation rules for the attribute.
        df.loc[df['Validation Rules'] == "", "Validation Rules"] = "_None_"

        # Drop the Valid Values column before final output
        df = df.drop(columns=['Valid Values'])

        df[COLS_TO_RENDER].to_csv(join(parent, 'template.csv'), index=False)