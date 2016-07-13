from pyramid import testing

import unittest
import mock
import json
import test.fixtures as factory
from test.mock_helpers import MockQuery
from src.views import search

class SearchTest(unittest.TestCase):    
    def setUp(self):
        testing.setUp()

        self.category_filters = {
        "locus": [('feature type', 'feature_type'), ('molecular function', 'molecular_function'), ('phenotype', 'phenotypes'), ('cellular component', 'cellular_component'), ('biological process', 'biological_process')],
        "phenotype": [("observable", "observable"), ("qualifier", "qualifier"), ("references", "references"), ("phenotype_locus", "phenotype_loci"), ("chemical", "chemical"), ("mutant_type", "mutant_type")],
        "biological_process": [("go_locus", "go_loci")],
        "cellular_component": [("go_locus", "go_loci")],
        "molecular_function": [("go_locus", "go_loci")],
        "reference": [("author", "author"), ("journal", "journal"), ("year", "year"), ("reference_locus", "reference_loci")]
    }

    
    def tearDown(self):
        testing.tearDown()
    
    @mock.patch('src.models.ESearch.search')
    def test_should_return_all_results_for_no_query_param_on_search(self, mock_es):
        request = testing.DummyRequest()
        request.context = testing.DummyResource()
        request.registry.settings['elasticsearch.index'] = 'searchable_items'
        response = search(request)

        mock_es.assert_called_with(body={
            'query': {'match_all': {}},
            'aggs': {
                'feature_type': {
                    'terms': {
                        'field': 'feature_type'
                    }
                },
                'categories': {
                    'terms': {
                        'field': 'category'
                    }
                }
            }, 'size': 0}, index='searchable_items')

    @mock.patch('src.models.ESearch.search')
    def test_should_apply_pagination(self, mock_es):
        request = testing.DummyRequest(params={'q': 'eu gene', 'limit': 10, 'offset': 25, 'category': 'unknown'})
        request.context = testing.DummyResource()
        request.registry.settings['elasticsearch.index'] = 'searchable_items'
        response = search(request)

        mock_es.assert_called_with(body={
            'highlight': {
                'fields': {
                    'observable': {}, 'mutant_type': {}, 'name': {}, 'author': {}, 'year': {}, 'journal': {}, 'molecular_function': {}, 'cellular_component': {}, 'phenotypes': {}, 'qualifier': {}, 'summary': {}, 'chemical': {}, 'go_loci': {}, 'references': {}, 'reference_loci': {}, 'phenotype_loci': {}, 'biological_process': {}, 'name_description': {}, 'description': {}}},
            'query': {
                'filtered': {
                    'filter': {
                        'bool': {
                            'must': [{'term': {'category': 'unknown'}}]
                        }
                    },
                    'query': {
                        'bool': {
                            'minimum_should_match': 1,
                            'should': [{
                                'match_phrase_prefix': {
                                    'name': {
                                        'query': 'eu gene',
                                        'boost': 4,
                                        'analyzer': 'standard',
                                        'max_expansions': 30
                                    }
                                }
                            }, {
                                'match_phrase': {
                                    'name': {
                                        'query': 'eu gene',
                                        'boost': 10,
                                        'analyzer': 'standard'
                                    }
                                }
                            }, {
                                'match': {
                                    'description': {
                                        'query': 'eu gene',
                                        'boost': 3,
                                        'analyzer': 'standard'
                                    }
                                }
                            }, {
                                'match_phrase': {
                                    'keys': {
                                        'query': 'eu gene',
                                        'boost': 20,
                                        'analyzer': 'standard'
                                    }
                                }
                            }, {
                                'multi_match': {
                                    'query': 'eu gene',
                                    'boost': 3,
                                    'type': 'best_fields',
                                    'fields': ['summary', 'name_description', 'phenotypes', 'cellular_component', 'biological_process', 'molecular_function', 'observable', 'qualifier', 'references', 'phenotype_loci', 'chemical', 'mutant_type', 'go_loci', 'author', 'journal', 'year', 'reference_loci']
                                }
                            }]
                        }
                    }
                }}, '_source': ['name', 'href', 'description', 'category', 'bioentity_id', 'keys']}, from_=25, index='searchable_items', size=10)

    @mock.patch('src.models.ESearch.search')
    @mock.patch('src.views.format_aggregation_results')
    @mock.patch('src.views.format_search_results')
    def test_should_format_response(self, mock_format_search, mock_format_aggs, mock_es):
        mock_es.return_value = {
            'hits': {
                'total': 10
            }
        }

        mock_format_search.return_value = [1,2,3]
        mock_format_aggs.return_value = [3,2,1]
        
        request = testing.DummyRequest(params={'limit': 10, 'offset': 25, 'category': 'locus'})
        request.context = testing.DummyResource()
        request.registry.settings['elasticsearch.index'] = 'searchable_items'
        response = search(request)

        self.assertEqual(response, {
            'total': 10,
            'results': [1,2,3],
            'aggregations': [3,2,1]
        })

    @mock.patch('src.models.ESearch.search')
    @mock.patch('src.views.format_aggregation_results')
    @mock.patch('src.views.format_search_results')
    def test_should_format_response(self, mock_format_search, mock_format_aggs, mock_es):
        mock_es.return_value = {
            'hits': {
                'total': 10
            }
        }

        mock_format_search.return_value = [1,2,3]
        mock_format_aggs.return_value = [3,2,1]
        
        request = testing.DummyRequest(params={'limit': 10, 'offset': 25, 'category': 'locus'})
        request.context = testing.DummyResource()
        request.registry.settings['elasticsearch.index'] = 'searchable_items'
        response = search(request)

        mock_format_aggs.assert_called_with(mock_es.return_value, 'locus', self.category_filters)
        
        self.assertEqual(response, {
            'total': 10,
            'results': [1,2,3],
            'aggregations': [3,2,1]
        })

    @mock.patch('src.models.ESearch.search')
    @mock.patch('src.views.format_aggregation_results')
    @mock.patch('src.views.format_search_results')
    def test_should_return_empty_aggregations_for_unknown_categories(self, mock_format_search, mock_format_aggs, mock_es):
        mock_es.return_value = {
            'hits': {
                'total': 10
            }
        }

        mock_format_search.return_value = [1,2,3]
        mock_format_aggs.return_value = [3,2,1]
        
        request = testing.DummyRequest(params={'limit': 10, 'offset': 25, 'category': 'blablabla'})
        request.context = testing.DummyResource()
        request.registry.settings['elasticsearch.index'] = 'searchable_items'
        response = search(request)

        mock_format_aggs.assert_called_with([], 'blablabla', self.category_filters)