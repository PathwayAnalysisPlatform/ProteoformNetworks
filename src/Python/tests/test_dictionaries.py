import os
from unittest import TestCase

from lib.dictionaries import read_dictionary_one_to_set, write_dictionary_one_to_one, convert_tab_to_dict, \
    write_dictionary_one_to_set, read_set_from_columns, get_intersection, in_dictionary, flatten_dictionary


class Test_dictionaries(TestCase):

    def setUp(self):
        self.letters = {chr(ord('A') + i): i for i in range(23)}
        languages = [('C++', '.cpp'), ('Python', '.py'), ('JavaScript', '.js'), ('C++', '.hpp'), ('C++', '.cpp')]
        with open('languages.txt', 'w') as file_languages:
            for entry in languages:
                file_languages.write(f"{entry[0]}\t{entry[1]}\n")
        self.languages = read_dictionary_one_to_set('./', 'languages.txt')

    def tearDown(self):
        os.remove('languages.txt')

    def test_write_dictionary_one_to_one_creates_file(self):
        letters = {chr(ord('A') + i): i for i in range(23)}
        file_name = 'test_write_dictionary_creates_file_alphabet.txt'
        file_path = './'
        self.assertFalse(os.path.exists(file_path + file_name))
        write_dictionary_one_to_one(self.letters, file_path, file_name)
        self.assertTrue(os.path.exists(file_path + file_name))
        os.remove(file_path + file_name)

    def test_write_dictionary_one_to_one_empty_path(self):
        file_name = 'test_write_dictionary_empty_path_alphabet.txt'
        file_path = ''
        self.assertFalse(os.path.exists(file_path + file_name))
        write_dictionary_one_to_one(self.letters, file_path, file_name)
        try:
            self.assertTrue(os.path.exists(file_path + file_name))
        finally:
            os.remove(file_path + file_name)

    def test_write_dictionary_one_to_one_path_does_not_exist(self):
        file_name = 'test_write_dictionary_path_does_not_exist_alphabet.txt'
        file_path = 'my/new/path/'
        self.assertFalse(os.path.exists(file_path + file_name))
        write_dictionary_one_to_one(self.letters, file_path, file_name)
        self.assertTrue(os.path.exists(file_path + file_name))
        os.remove(file_path + file_name)
        os.rmdir('my/new/path/')
        os.rmdir('my/new/')
        os.rmdir('my/')

    def test_read_dictionary_reads_all_keys(self):
        self.assertIn('C++', self.languages.keys())
        self.assertIn('Python', self.languages.keys())
        self.assertIn('JavaScript', self.languages.keys())

    def test_read_dictionary_values_are_sets(self):
        self.assertEqual(set, type(self.languages['C++']))
        self.assertEqual(set, type(self.languages['Python']))

    def test_read_dictionary_removes_duplicate_entries(self):
        self.assertEqual(3, len(self.languages.keys()), msg=f'There is a repeated key: {self.languages.keys()}')
        self.assertEqual(2, len(self.languages['C++']),
                         msg=f"There is a duplicated value for the 'C++' key: {self.languages['C++']}")

    def test_read_dictionary_order_pairs_true(self):
        # Create file with pairs. Some with inverted lexicographic order
        pairs = [('a', 'b'), ('c', 'b'), ('d', 'e')]
        file_name = TestCase.id(self) + '_pairs.txt'
        with open(file_name, 'w') as file:
            for x, y in pairs:
                file.write(f"{x}\t{y}\n")
        # Execute target method
        result = read_dictionary_one_to_set('', file_name, order_pairs=True)

        # Check the pairs order was corrected, showing them as key and value when word1 < word2 Lexicographical order
        self.assertIn('b', result.keys(), msg="Missing key because it did not order the column values")
        self.assertEqual(3, len(result.keys()), msg="Wrong number of columns")

        os.remove(file_name)

    def test_read_dictionary_indices_1_2(self):
        # Create file with three columns, some not in lexicographic order
        trios = [(1, 1, 2), (2, 3, 2), (3, 4, 5)]
        file_name = TestCase.id(self) + '_pairs.txt'
        with open(file_name, 'w') as file:
            for x, y, z in trios:
                file.write(f"{x}\t{y}\t{z}\n")

        # Execute target method
        result = read_dictionary_one_to_set('', file_name, order_pairs=True, col_indices=(1, 2))

        # Check values are correct
        self.assertIn('1', result.keys(), msg="Missing key in dictionary")
        self.assertIn('2', result.keys(), msg="Missing key in dictionary")
        self.assertNotIn('3', result.keys(), msg="Incorrect key in dictionary")
        self.assertIn('4', result.keys(), msg="Missing key in dictionary")

        # Remove file
        os.remove(file_name)

    def test_read_dictionary_missing_two_columns(self):
        """With a one column file, request default columns 0 and 1, report error"""
        # Create file with three columns, some not in lexicographic order
        file_name = TestCase.id(self) + '_single_column.txt'
        with open(file_name, 'w') as file:
            for x in range(5):
                file.write(f"{x}\n")

        with self.assertRaises(ValueError,
                               msg='Should raise an exception because needed columns of the file are missing.'):
            read_dictionary_one_to_set('', file_name, order_pairs=True)

        os.remove(file_name)

    def test_read_dictionary_missing_index_columns(self):
        """With two columns file, indices other than (0, 1), like (1, 2), show error."""
        # Create file with three columns, some not in lexicographic order
        pairs = [('a', 'b'), ('c', 'b'), ('d', 'e')]
        file_name = TestCase.id(self) + '_pairs.txt'
        with open(file_name, 'w') as file:
            for x, y in pairs:
                file.write(f"{x}\t{y}\n")

        with self.assertRaises(ValueError,
                               msg='Should raise an exception because needed columns of the file are missing.'):
            read_dictionary_one_to_set('', file_name, col_indices=(1, 2))

        os.remove(file_name)

    def test_read_dictionary_skip_header(self):
        # Create trio file with headers
        trios = [('Column1', 'Column2', 'Column3'), (1, 1, 2), (2, 3, 2), (3, 4, 5)]
        file_name = TestCase.id(self) + '_pairs.txt'
        with open(file_name, 'w') as file:
            for x, y, z in trios:
                file.write(f"{x}\t{y}\t{z}\n")

        # Execute target method
        result = read_dictionary_one_to_set('', file_name, order_pairs=True, col_indices=(1, 2), ignore_header=True)

        # Check headers are not taken as key, value pairs
        self.assertNotIn('Column1', result.keys(), msg="Missing key in dictionary")
        self.assertIn('1', result.keys(), msg="Missing key in dictionary")
        self.assertIn('2', result.keys(), msg="Missing key in dictionary")
        self.assertIn('4', result.keys(), msg="Missing key in dictionary")

        # Remove precondition files
        os.remove(file_name)


class Test_convert_tab_to_dict(TestCase):

    def setUp(self):
        self.response = f"language\tending\nC++\t.cpp\nPython\t.py\nJavaScript\t.js\nC++\t.hpp\nC++\t.cpp"
        self.languages = convert_tab_to_dict(self.response)

    def test_convert_tab_to_dict_values_are_sets(self):
        self.assertEqual(set, type(self.languages['C++']), msg='The values are not sets')
        self.assertEqual(set, type(self.languages['Python']), msg='The values are not sets')

    def test_convert_tab_to_dict_correct_keys(self):
        self.assertEqual(3, len(self.languages.keys()), msg='Wrong number of keys')
        self.assertIn('C++', self.languages.keys(), msg='Missing dictionary key')
        self.assertIn('Python', self.languages.keys(), msg='Missing dictionary key')
        self.assertIn('JavaScript', self.languages.keys(), msg='Missing dictionary key')

    def test_convert_tab_to_dict_unique_values(self):
        self.assertEqual(2, len(self.languages['C++']),
                         msg=f"Wrong number of values for the key C++: {self.languages['C++']}")

    def test_convert_tab_to_dict_ignores_header(self):
        self.assertNotIn('language', self.languages.keys(), msg='Should not include the header line of the response')

    def test_convert_tab_to_dict_not_string_returns_empty_dict(self):
        self.assertEqual({}, convert_tab_to_dict(5.5), msg='Should return an empty dictionary')


class Test_write_dictionary_one_to_set(TestCase):

    def setUp(self):
        self.friends = {'David': {'Ester', 'Martha', 'Mathew'},
                        'Michael': {'David'},
                        'James': {'Maria', 'Martha'}}

    def test_write_dictionary_one_to_set_non_existent_path_creates_file(self):
        file_name = TestCase.id(self) + '_friends.txt'
        file_path = 'a/path/'
        self.assertFalse(os.path.exists(file_path + file_name), msg='The file should not exist as precondition')
        write_dictionary_one_to_set(self.friends, file_path, file_name)
        self.assertTrue(os.path.exists(file_path + file_name), msg='The file was not created')
        os.remove(file_path + file_name)
        os.rmdir('a/path/')
        os.rmdir('a/')

    def test_write_dictionary_one_to_set_empty_path_works(self):
        file_name = TestCase.id(self) + '_friends.txt'
        self.assertFalse(os.path.exists(file_name), msg='The file should not exist as precondition')
        write_dictionary_one_to_set(self.friends, '', file_name)
        self.assertTrue(os.path.exists(file_name), msg='The file was not created')
        os.remove(file_name)

    def test_write_dictionary_one_to_set_empty_file_name_raises_error(self):
        with self.assertRaises(ValueError, msg="Should not accept empty file name"):
            write_dictionary_one_to_set({}, '', '')

    def test_write_dictionary_one_to_set_one_row_each_set_entry(self):
        """Having a key with multiple values, writes one row for each key -- value pair"""
        file_name = TestCase.id(self) + '_friends.txt'
        self.assertFalse(os.path.exists(file_name), msg='The file should not exist as precondition')
        write_dictionary_one_to_set(self.friends, '', file_name)

        with open(file_name) as file:
            lines = file.readlines()
            self.assertEqual(6, len(lines), msg='Wrong number of lines in the file')
            self.assertIn('David\tEster\n', lines, msg='Missing entry line')
            self.assertIn('David\tMartha\n', lines, msg='Missing entry line')
            self.assertIn('David\tMathew\n', lines, msg='Missing entry line')
        os.remove(file_name)

    def test_write_dictionary_one_to_set_all_keys(self):
        file_name = TestCase.id(self) + '_friends.txt'
        self.assertFalse(os.path.exists(file_name), msg='The file should not exist as precondition')
        write_dictionary_one_to_set(self.friends, '', file_name)

        with open(file_name) as file:
            lines = file.readlines()
            self.assertIn('David\tEster\n', lines, msg='Missing entry line')
            self.assertIn('Michael\tDavid\n', lines, msg='Missing entry line')
            self.assertIn('James\tMartha\n', lines, msg='Missing entry line')
        os.remove(file_name)


class Test_read_set_column_values(TestCase):

    def setUp(self):
        # Create file with trios, no header
        trios = [(1, 'A', 'B'), (2, 'B', 'C'), (3, 'D', 'E'), (3, 'D', 'E'), ('B', 2, 'E')]
        self.file_name_trios = 'trios.txt'
        with open(self.file_name_trios, 'w') as file:
            for x, y, z in trios:
                file.write(f"{x}\t{y}\t{z}\n")
        self.result_trios = read_set_from_columns('', self.file_name_trios)

    def tearDown(self):
        os.remove(self.file_name_trios)

    def test_read_set_columns_correct_result_type(self):
        self.assertEqual(set, type(self.result_trios), msg='Wrong return type of the function')

    def test_read_set_columns_removes_duplicates(self):
        # Removes duplicates in the same column and mixing
        self.assertEqual(6, len(self.result_trios), msg='Wrong number of values in set')
        unique_values = ['1', '2', '3', 'A', 'B', 'D']  # From columns 0 and 1
        for value in unique_values:
            self.assertIn(value, self.result_trios, msg='Missing value in the set')

    def test_read_set_from_columns_ignores_header(self):
        result_ignore_header = read_set_from_columns('', self.file_name_trios, ignore_header=True)
        self.assertEqual(4, len(result_ignore_header), msg='Wrong number of values in set')
        unique_values = ['2', '3', 'B', 'D']  # From columns 0 and 1
        for value in unique_values:
            self.assertIn(value, result_ignore_header, msg='Missing value in the set')

    def test_read_set_from_columns_non_existent_indices_raises_exception(self):
        with self.assertRaises(ValueError, msg='Should raise exception for inexistent column'):
            read_set_from_columns('', self.file_name_trios, col_indices=(0, 1, 2, 3))

    def test_read_set_from_columns_no_columns_raises_exception(self):
        with self.assertRaises(ValueError, msg='Should raise exception for missing columns argument'):
            read_set_from_columns('', self.file_name_trios, col_indices=())

    def test_read_set_from_columns_one_column_selected(self):
        result_one_column = read_set_from_columns('', self.file_name_trios, col_indices=[0])
        self.assertEqual(4, len(result_one_column), msg='Wrong number of values')
        self.assertIn('B', result_one_column, msg='Missing value in the set')
        self.assertNotIn('A', result_one_column, msg='The key of another column should not be added')

    def test_read_set_from_columns_multiple_columns_selected(self):
        result_three_columns = read_set_from_columns('', self.file_name_trios, col_indices=[0, 1, 2])
        self.assertEqual(8, len(result_three_columns), msg=f"Wrong number of values: {result_three_columns}")
        self.assertIn('B', result_three_columns, msg='Missing value in the set')
        self.assertIn('A', result_three_columns, msg='Missing value in the set')
        self.assertIn('E', result_three_columns, msg='Missing value from column in the set')

    def test_get_intersection(self):
        # Create dictionary1
        dictionary1 = {'A': {'a', 'b', 'c'},
                       'B': {'b'},
                       'C': {'c', 'd', 'e'}}
        # Create dictionary2
        dictionary2 = {'A': {'a', 'c'},
                       'B': {'c'},
                       'C': {'c'}}

        intersection = get_intersection(dictionary1, dictionary2)
        # Get pairs of dictionary1 that appear in dictionary2
        self.assertIn('A', intersection.keys(), msg='Missing key')
        self.assertNotIn('B', intersection.keys(), msg='Key should not be in the dictionary')
        self.assertIn('C', intersection.keys(), msg='Missing key')
        self.assertIn('a', intersection['A'], msg=f"Missing value for the key 'A'")
        self.assertIn('c', intersection['A'], msg=f"Missing value for the key 'A'")
        self.assertIn('c', intersection['C'], msg=f"Missing value for the key 'C'")
        self.assertNotIn('d', intersection['C'], msg=f"Key should not be a value of 'C'")

    def test_flatten_dictionary(self):
        dictionary = {'A': {'a', 'b', 'c'},
                       'B': {'b'},
                       'C': {'c', 'd', 'e'}}
        result = flatten_dictionary(dictionary)
        self.assertEqual(7, len(result), msg='Wrong number of elements in the list')
        self.assertIn(('A', 'a'), result, msg='Missing entry in the list')
        self.assertIn(('A', 'c'), result, msg='Missing entry in the list')
        self.assertIn(('B', 'b'), result, msg='Missing entry in the list')
        self.assertIn(('C', 'c'), result, msg='Missing entry in the list')
        self.assertIn(('C', 'e'), result, msg='Missing entry in the list')

    def test_in_dictionary(self):
        # Create list of pairs
        pairs = [('A', 'B'), ('B', 'C'), ('D', 'C'), ('E', 'D')]
        # Create dictionary
        dictionary = {'A': {'B', 'C'},
                      'D': {'C'},
                      'E': {'D'}}

        # Return a boolean vector saying if each pair is contained in the dictionary
        result = in_dictionary(pairs, dictionary)

        # Assert correctness
        self.assertEqual(4, len(result), msg=f"The result list should be {len(pairs)} elements long")
        self.assertTrue(result[0], msg="The pair should be marked as contained in the dictionary")
        self.assertFalse(result[1], msg="The pair should not be marked as contained in the dictionary")
        self.assertTrue(result[2], msg="The pair should be marked as contained in the dictionary")
        self.assertTrue(result[3], msg="The pair should be marked as contained in the dictionary")