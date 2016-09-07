def add_sort_by(sort_by, search_body):
    if sort_by == 'alphabetical':
        search_body['sort'] = [
            {
                "name.raw": {
                    "order": "asc"
                }
            }
        ]
    elif sort_by == 'annotation':
        search_body['sort'] = [
            {
                "number_annotations": {
                    "order": "desc"
                }
            }
        ]

def add_highlighting(fields, search_body):
    for field in fields:
        search_body['highlight']['fields'][field] = {}

def format_search_results(search_results, query, response_fields):
    formatted_results = []

    for r in search_results['hits']['hits']:
        raw_obj = r.get('_source')

        obj = {}
        for field in response_fields:
            obj[field] = raw_obj.get(field)
                
        obj['highlights'] = r.get('highlight')

        formatted_results.append(obj)

    if search_is_quick(query, search_results):
        formatted_results[0]['is_quick'] = True

    return formatted_results

def format_aggregation_results(aggregation_results, category, category_filters):
    if category == '':
        category_obj = {'values': [], 'key': 'category'}
        for bucket in aggregation_results['aggregations']['categories']['buckets']:
            category_obj['values'].append({'key': bucket['key'], 'total': bucket['doc_count']})
        return [category_obj]
    
    elif category in category_filters.keys():
        formatted_agg = []
        
        for agg_info in category_filters[category]:
            agg_obj = {'key': agg_info[0], 'values': []}
            if agg_info[1] in aggregation_results['aggregations']:
                for agg in aggregation_results['aggregations'][agg_info[1]]['buckets']:
                    agg_obj['values'].append({'key': agg['key'], 'total': agg['doc_count']})
            formatted_agg.append(agg_obj)
            
        return formatted_agg
    
    else:
        return []

def search_is_quick(query, search_results):
    if search_results['hits']['hits'][0].get('_source').get('keys'):
        if query and query.replace('"', '').lower().strip() in search_results['hits']['hits'][0].get('_source').get('keys'):
            if len(search_results['hits']['hits']) > 1:
                if (query.replace('"', '').lower().strip() not in search_results['hits']['hits'][1].get('_source').get('keys')):
                    return True
            else:
                return True
    return False

def add_exact_match(es_query, query, multi_match_fields):
    new_conditions = []
    for cond in es_query['bool']['should'][2:4]:
        new_conditions.append({'match_phrase_prefix': cond.pop(cond.keys()[0])})
    multi_fields = {
        "multi_match": {
            "query": query,
            "type": "phrase_prefix",
            "fields": multi_match_fields,
            "boost": 3
        }
    }
    es_query['bool']['should'] = [es_query['bool']['should'][0]] + new_conditions + [multi_fields]

def build_aggregation_query(es_query, category, category_filters):
    agg_query_body = {
        'query': es_query,
        'size': 0,
        'aggs': {}
    }

    if category == '':
        agg_query_body['aggs'] = {
            'categories': {
                'terms': {'field': 'category', 'size': 50}
            },
            'feature_type': {
                'terms': {'field': 'feature_type', 'size': 50}
            }
        }
    else:
        for c in category_filters[category]:
            agg_query_body['aggs'][c[1]] = {'terms': {'field': c[1] + '.raw', 'size': 999}}

    return agg_query_body

def build_es_search_body(query, category, es_query, returning_fields):
    results_search_body = {
        '_source': returning_fields,
        'highlight': {
            'fields' : {}
        },
        'query': es_query
    }
    
    if query == '' and category == '':
        results_search_body["query"] = {
            "function_score": {
                "query": es_query,
                "random_score": {"seed" : 12345}
            }
        }
    else:
        results_search_body["sort"] = [
            '_score',
            {'number_annotations': {'order': 'desc'}}
        ]

    return results_search_body

def build_search_query(query, multi_match_fields, category, category_filters, request):
    es_query = build_es_search_query(query, multi_match_fields)

    if category == '':
        return es_query
    
    query = {
        'filtered': {
            'query': es_query,
            'filter': {
                'bool': {
                    'must': [{'term': {'category': category}}]
                }
            }
        }
    }
    
    if category in category_filters.keys():
        for item in category_filters[category]:
            if request.params.get(item[0]):
                query['filtered']['filter']['bool']['must'].append({'term': {(item[1]+".raw"): request.params.get(item[0])}})

    return query
    
def build_es_search_query(query, multi_match_fields):
    for special_char in ['-', '.']:
        if special_char in query:
            query = "\"" + query + "\""
            break
    
    if query == '':
        return {'match_all': {}}
    
    es_query = {
        "bool": {
            "should": [
                {
                    "match_phrase_prefix": {
                        "name": {
                            "query": query,
                            "boost": 3,
                            "max_expansions": 30,
                            "analyzer": "standard"
                        }
                    }
                },
                {
                    "match_phrase_prefix": {
                        "keys": {
                            "query": query,
                            "boost": 35,
                            "max_expansions": 12,
                            "analyzer": "standard"
                        }
                    }
                },
                {
                    "match_phrase": {
                        "name": {
                            "query": query,
                            "boost": 80,
                            "analyzer": "standard"
                        }
                    }
                },                        
                {
                    "match": {
                        "description": {
                            "query": query,
                            "boost": 1,
                            "analyzer": "standard"
                        }
                    }
                },
                {
                    "match_phrase": {
                        "keys": {
                            "query": query,
                            "boost": 50,
                            "analyzer": "standard"
                        }
                    }
                },
                {
                    "multi_match": {
                        "query": query,
                        "type": "best_fields",
                        "fields": multi_match_fields,
                        "boost": 1
                    }
                }
            ],
            "minimum_should_match": 1
        }
    }

    if (query[0] in ('"', "'") and query[-1] in ('"', "'")):
        add_exact_match(es_query, query, multi_match_fields)

    return es_query

def build_autocomplete_search_query(query, category):
    es_query = {
        "query": {
            "bool": {
                "must": [{
                    "match": {
                        "name": {
                            "query": query,
                            "analyzer": "standard"
                        }
                    }
                }],
                "must_not": [{
                    "match": { "category": "reference" }
                }, {
                    "match": { "category": "download" }
                }],
                "should": [
                    {
                        "match": {
                            "category": {
                                "query": "locus",
                                "boost": 4
                            }
                        }
                    }
                ]
            }
        }, '_source': ['name', 'href', 'category']
    }

    if category != '':
        es_query["query"]["bool"]["must"].append({"match": {"category": category}})
        if category != "locus":
            es_query["query"]["bool"].pop("should")

    return es_query

def format_autocomplete_results(es_response):
    formatted_results = []
    for hit in es_response['hits']['hits']:
        obj = {
            'name': hit['_source']['name'],
            'href': hit['_source']['href'],
            'category': hit['_source']['category']
        }
        formatted_results.append(obj)
    return formatted_results
