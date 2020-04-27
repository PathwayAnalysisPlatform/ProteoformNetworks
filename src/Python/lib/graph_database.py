from neo4j import GraphDatabase
import pandas as pd

def get_query_result(query):
    """
    Get pandas dataframe with the result of the query to Reactome Graph Database
    :param query: Cypher string with the query to get the data
    :return: pandas dataframe with the result records
    """
    db = GraphDatabaseAccess("bolt://localhost", "neo4j", "")
    df = db.get_result(query)
    db.close()
    return df;

class GraphDatabaseAccess:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password), encrypted=False)

    def close(self):
        self.driver.close()

    def get_result(self, query):
        with self.driver.session() as session:
            records = session.read_transaction(self.get_records, query)
            df = pd.DataFrame([r.values() for r in records], columns=records[0].keys())
            return df

    @staticmethod
    def get_records(tx, query):
        result = []
        for record in tx.run(query):
            result.append(record)
        return result

if __name__ == "__main__":
    pathway = "R-HSA-70171"
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId"
    df = get_query_result(query)
    print(df)
