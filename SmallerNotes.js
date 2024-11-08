//Smaller things:
//You can use Infinity as an actual pre-declared value. 


//

// //1
// function make_k_list(k, d) {
// if (d === 0) {
// return 0;
// } else {
// let klist = null;
// for (let i = 0; i < k; i = i + 1) {
// klist = pair(make_k_list(k, d - 1), klist);
// }
// return klist;
// }
// }
// //===============================================================
// // TASK 1B
// //===============================================================
// function sum_k_list(klist) {
// if (is_number(klist)) {
// return klist;
// } else {
// const k = length(klist);
// let sum = 0;
// let rest = klist;
// for (let i = 0; i < k; i = i + 1) {
// sum = sum + sum_k_list(head(rest));
// rest = tail(rest);
// }
// return sum;
// }
// }
// //===============================================================
// // TASK 1C
// //===============================================================
// function map_k_list(f, klist) {
// if (is_number(klist)) {
// return f(klist);
// } else {
// const k = length(klist);
// let result = null;
// let rest = klist;
// for (let i = 0; i < k; i = i + 1) {
// result = pair(map_k_list(f, head(rest)), result);
// rest = tail(rest);
// }
// return reverse(result);
// }
// }
// //===============================================================
// // TASK 2A
// //===============================================================
// function route_distance(mat, route) {
// function add_dist(rou) {
// if (is_null(rou) || is_null(tail(rou))) {
// return 0;
// } else {
// return mat[head(rou)][head(tail(rou))] + add_dist(tail(rou));
// }
// }
// return add_dist(route);
// }
// //===============================================================
// // TASK 2B
// //===============================================================
// function shortest_paper_route(n, mat, start) {
// // You can keep, modify or remove the permutations function.
// function permutations(ys) {
// return is_null(ys)
// ? list(null)
// : accumulate(append, null,
// map(x => map(p => pair(x, p),
// permutations(remove(x, ys))),
// ys));
// }
// const others = remove(start, enum_list(0, n - 1));
// const routes = permutations(others);
// let min_dist = Infinity;
// let min_route = null;
// for (let p = routes; !is_null(p); p = tail(p)) {
// const pp_route = pair(start, append(head(p), list(start)));
// const route_dist = route_distance(mat, pp_route);
// if (route_dist < min_dist) {
// min_dist = route_dist;
// min_route = pp_route;
// } else { }
// }
// return pair(min_route, min_dist);
// }
// //===============================================================
// // TASK 3A
// //===============================================================
// function make_postfix_exp(bae) {
// const pfe = [];
// let next = 0;
// function convert(sub_bae) {
// if (is_number(sub_bae)) {
// pfe[next] = sub_bae;
// next = next + 1;
// } else {
// convert(sub_bae[0]);
// convert(sub_bae[2]);
// pfe[next] = sub_bae[1];
// next = next + 1;
// }
// }
// convert(bae);
// return pfe;
// }
// //===============================================================
// // TASK 3B
// //===============================================================
// function eval_postfix_exp(pfe) {
// let next = array_length(pfe) - 1;
// function evaluate() {
// const token = pfe[next];
// next = next - 1;
// if (is_number(token)) {
// return token;
// } else {
// const op = token;
// const right = evaluate();
// const left = evaluate();
// if (op === "+") {
// return left + right;
// } else if (op === "-") {
// return left - right;
// } else if (op === "*") {
// return left * right;
// } else {
// return left / right;
// }
// }
// }
// return evaluate();
// }





// //===============================================================
// // Task 1A: Delta Encoding
// //===============================================================
// function delta_encode(L) {
// function encode(xs, prev) {
// return is_null(xs)
// ? null
// : pair(head(xs) - prev, encode(tail(xs), head(xs)));
// }
// return encode(L, 0);
// }
// //===============================================================
// // Task 1B: Delta Decoding
// //===============================================================
// function delta_decode(D) {
// function decode(xs, prev) {
// return is_null(xs)
// ? null
// : pair(prev + head(xs), decode(tail(xs), prev + head(xs)));
// }
// return decode(D, 0);
// }
// //===============================================================
// // Task 2A: Run-Length Encoding
// //===============================================================
// function runlength_encode(L) {
// function encode(val, count, next) {
// return is_null(next)
// ? list(count === 1 ? val : pair(val, count))
// : head(next) === val
// ? encode(val, count + 1, tail(next))
// : pair(count === 1 ? val : pair(val, count),
// encode(head(next), 1, tail(next)));
// }
// return is_null(L)
// ? null
// : encode(head(L), 1, tail(L));
// }
// //===============================================================
// // Task 2B: Run-Length Decoding
// //===============================================================
// function runlength_decode(R) {
// function decode(xs, val, count) {
// return count > 0
// ? pair(val, decode(xs, val, count - 1))
// : is_null(xs)
// ? null
// : !is_pair(head(xs))
// ? pair(head(xs), decode(tail(xs), 0, 0))
// : decode(tail(xs), head(head(xs)), tail(head(xs)));
// }
// return decode(R, 0, 0);
// }
// //===============================================================
// // Task 3A: Smallest Bounding Rectangle
// //===============================================================
// const get_x = (aar) => list_ref(aar, 0);
// const get_y = (aar) => list_ref(aar, 1);
// const get_width = (aar) => list_ref(aar, 2);
// const get_height = (aar) => list_ref(aar, 3);
// function smallest_bounding_AAR_area(rs) {
// let min_x = Infinity;
// let min_y = Infinity;
// let max_x = -Infinity;
// let max_y = -Infinity;
// for (let p = rs; !is_null(p); p = tail(p)) {
// const aar = head(p);
// const x1 = get_x(aar);
// const x2 = x1 + get_width(aar);
// const y1 = get_y(aar);
// const y2 = y1 + get_height(aar);
// if (x1 < min_x) { min_x = x1; } else { }
// if (x2 > max_x) { max_x = x2; } else { }
// if (y1 < min_y) { min_y = y1; } else { }
// if (y2 > max_y) { max_y = y2; } else { }
// }
// return (max_x - min_x) * (max_y - min_y);
// }
// //===============================================================
// // Task 3B: Optimized Smallest Bounding Rectangle
// //===============================================================
// const get_x = (aar) => list_ref(aar, 0);
// const get_y = (aar) => list_ref(aar, 1);
// const get_width = (aar) => list_ref(aar, 2);
// const get_height = (aar) => list_ref(aar, 3);
// function optimized_smallest_bounding_AAR_area(rs) {
// let max_longer = 0;
// let max_shorter = 0;
// for (let p = rs; !is_null(p); p = tail(p)) {
// const aar = head(p);
// const width = get_width(aar);
// const height = get_height(aar);
// const longer = math_max(width, height);
// const shorter = math_min(width, height);
// if (longer > max_longer) { max_longer = longer; } else { }
// if (shorter > max_shorter) { max_shorter = shorter; } else { }
// }
// return max_longer * max_shorter;
// }
// //===============================================================
// // Task 3C: Overlapping Rectangles
// //===============================================================
// const get_x = (aar) => list_ref(aar, 0);
// const get_y = (aar) => list_ref(aar, 1);
// const get_width = (aar) => list_ref(aar, 2);
// const get_height = (aar) => list_ref(aar, 3);
// // SOLUTION 1
// function overlap_area(aar1, aar2) {
// // [a, b] and [c, d] are the input intervals.
// function overlap_length(a, b, c, d) {
// return math_max(0, math_min(b, d) - math_max(a, c));
// }
// const x_overlap = overlap_length(
// get_x(aar1), get_x(aar1) + get_width(aar1),
// get_x(aar2), get_x(aar2) + get_width(aar2));
// const y_overlap = overlap_length(
// get_y(aar1), get_y(aar1) + get_height(aar1),
// get_y(aar2), get_y(aar2) + get_height(aar2));
// return x_overlap * y_overlap;
// }
// // SOLUTION 2
// function overlap_area(aar1, aar2) {
// // [a, b] and [c, d] are the input intervals.
// function overlap_length(a, b, c, d) {
// if (c < a) {
// return overlap_length(c, d, a, b); // to make sure a <= c
// } else if (b <= c) {
// return 0; // when a < b <= c < d
// } else if (b <= d) {
// return b - c; // when a <= c < b <= d
// } else {
// return d - c; // when a <= c < d <= b
// }
// }
// const x_overlap = overlap_length(
// get_x(aar1), get_x(aar1) + get_width(aar1),
// get_x(aar2), get_x(aar2) + get_width(aar2));
// const y_overlap = overlap_length(
// get_y(aar1), get_y(aar1) + get_height(aar1),
// get_y(aar2), get_y(aar2) + get_height(aar2));
// return x_overlap * y_overlap;
// }




// //=====================================================================
// // TASK 1A
// function is_pa_word(s) {
// return !is_null(member(s, pa_words));
// }
// // Alternative solution
// function is_pa_word(s) {
// return !is_null(filter(x => x === s, pa_words));
// }
// // Alternative solution
// function is_pa_word(s) {
// return accumulate((x, y) => (x === s) || y, false, pa_words);
// }
// // Alternative solution
// function is_pa_word(s) {
// function iter(xs) {
// return is_null(xs) ? false : s === head(xs) ? true : iter(tail(xs));
// }
// return iter(pa_words);
// }
// //=====================================================================
// // TASK 1B
// function count_matches(char, pos) {
// return length(filter(x => char_at(x, pos) === char, pa_words));
// }
// //=====================================================================
// // TASK 1C
// function helper(s, i) {
// const c = char_at(s, i);
// return is_undefined(c)
// ? null
// : pair(c, () => helper(s, i + 1));
// }
// function char_stream(s) {
// return helper(s, 0);
// }
// //=====================================================================
// // TASK 1D
// function solve(n, constraints) {
// return accumulate((constraint, ss) =>
// filter(s => char_at(s, head(constraint))
// === tail(constraint), ss),
// filter(s => string_length(s) === n, pa_words),
// constraints);
// }
// // Alternative solution (less efficient)
// function solve(n, constraints) {
// function good_word(w, constraints) {
// return accumulate((c, good) => char_at(w, head(c)) === tail(c) ? good :
// false,
// true,
// constraints);
// }
// const pa_words_len_n = filter(s => string_length(s) === n, pa_words);
// return filter(w => good_word(w, constraints), pa_words_len_n);
// }
// //=====================================================================
// //=====================================================================
// // TASK 2A
// function eval_poly(poly) {
// function p(x) {
// return accumulate(
// (t, sum) => head(t) * math_pow(x, tail(t)) + sum,
// 0,
// poly);
// }
// return p;
// }
// // Alternative solution (using Horner's scheme)
// function eval_poly(poly) {
// function horner(poly, x, exp, acc) {
// return exp < 0
// ? acc
// : tail(head(poly)) < exp
// ? horner(poly, x, exp - 1, acc * x)
// : horner(tail(poly), x, exp - 1, head(head(poly)) + acc * x);
// }
// if (is_null(poly)) {
// return x => 0;
// } else {
// const rev_poly = reverse(poly);
// return x => horner(rev_poly, x, tail(head(rev_poly)), 0);
// }
// }
// //=====================================================================
// // TASK 2B
// function add_poly(poly1, poly2) {
// if (is_null(poly1)) {
// return poly2;
// } else if (is_null(poly2)) {
// return poly1;
// } else {
// const coeff1 = head(head(poly1));
// const coeff2 = head(head(poly2));
// const exp1 = tail(head(poly1));
// const exp2 = tail(head(poly2));
// if (exp1 === exp2) {
// return coeff1 + coeff2 === 0
// ? add_poly(tail(poly1), tail(poly2))
// : pair(pair(coeff1 + coeff2, exp1),
// add_poly(tail(poly1), tail(poly2)));
// } else if (exp1 < exp2) {
// return pair(head(poly1), add_poly(tail(poly1), poly2));
// } else {
// return pair(head(poly2), add_poly(poly1, tail(poly2)));
// }
// }
// }
// //=====================================================================
// // TASK 2C
// function multiply_poly(poly1, poly2) {
// const coeff = head;
// const exp = tail;
// return accumulate((term1, acc) =>
// add_poly(map(term2 =>
// pair(coeff(term1) * coeff(term2),
// exp(term1) + exp(term2)),
// poly2),
// acc),
// null,
// poly1);
// }
// // Alternative solution
// function multiply_poly(poly1, poly2) {
// return accumulate((p, q) => add_poly(p, q),
// null,
// map(t1 => map(t2 => pair(head(t1) * head(t2),
// tail(t1) + tail(t2)),
// poly2),
// poly1));
// }
// //=====================================================================
// //=====================================================================
// // TASK 3
// function alt_column_matrix(R, C) {
// const M = [];
// for (let r = 0; r < R; r = r + 1) {
// M[r] = [];
// }
// let count = 1;
// for (let c = 0; c < C; c = c + 1) {
// if (c % 2 === 0) {
// for (let r = 0; r < R; r = r + 1) {
// M[r][c] = count;
// count = count + 1;
// }
// } else {
// for (let r = R - 1; r >= 0; r = r - 1) {
// M[r][c] = count;
// count = count + 1;
// }
// }
// }
// return M;
// }
// //=====================================================================



// //=====================================================================
// // TASK 1A
// function split(S) {
// let i = 0;
// let result = [];
// while (char_at(S, i) !== undefined) {
// result[i] = char_at(S, i);
// i = i + 1;
// }
// return result;
// }
// //=====================================================================
// // TASK 1B
// function contains(B, A_i) {
// const lenB = array_length(B);
// for (let j = 0; j < lenB; j = j + 1) {
// if (B[j] === A_i) {
// return true;
// }
// }
// return false;
// }
// function num_characters_from(A, B) {
// const len = array_length(A);
// let count = 0;
// for (let i = 0; i < len; i = i + 1) {
// if (contains(B, A[i])) {
// count = count + 1;
// }
// }
// return count;
// }
// //=====================================================================
// // TASK 1C
// function num_unique(A) {
// const len = array_length(A);
// let count = 0;
// for (let i = 0; i < len; i = i + 1) {
// let last_occurrence = true;
// for (let j = i + 1; j < len; j = j + 1) {
// if (A[j] === A[i]) {
// last_occurrence = false;
// }
// }
// if (last_occurrence) {
// count = count + 1;
// }
// }
// return count;
// }
// //=====================================================================
// // TASK 2A
// function search_array(A, x) {
// const len = array_length(A);
// let i = 0;
// while (i < len && A[i] !== x) {
// i = i + 1;
// }
// return i < len ? i : -1;
// }
// function baseN_to_value(X) {
// const DIGIT_SYMBOLS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
// "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
// "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
// "U", "V", "W", "X", "Y", "Z"];
// const base = head(X);
// const digits = tail(X);
// let value = 0;
// for (let p = digits; !is_null(p); p = tail(p)) {
// value = value * base + search_array(DIGIT_SYMBOLS, head(p));
// }
// return value;
// }
// //=====================================================================
// // TASK 2B
// function value_to_baseN(N, x) {
// const DIGIT_SYMBOLS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
// "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
// "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
// "U", "V", "W", "X", "Y", "Z"];
// let remainder = x % N;
// let digits = list(DIGIT_SYMBOLS[remainder]);
// let quotient = math_floor(x / N);
// while (quotient > 0) {
// remainder = quotient % N;
// digits = pair(DIGIT_SYMBOLS[remainder], digits);
// quotient = math_floor(quotient / N);
// }
// return pair(N, digits);
// }
// // Alternative solution:
// function value_to_baseN(N, x) {
// const DIGIT_SYMBOLS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
// "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
// "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
// "U", "V", "W", "X", "Y", "Z"];
// function iter(val, digits) {
// const remainder = DIGIT_SYMBOLS[val % N];
// const quotient = math_floor(val / N);
// return quotient === 0
// ? pair(remainder, digits)
// : iter(quotient, pair(remainder, digits));
// }
// return pair(N, iter(x, null));
// }
// //=====================================================================
// // TASK 3A
// function flatten_bin_tree(T) {
// return is_null(T)
// ? null
// : append(flatten_bin_tree(list_ref(T, 1)),
// pair(head(T), flatten_bin_tree(list_ref(T, 2))));
// }
// //=====================================================================
// // TASK 3B
// function insert(x, xs) {
// return is_null(xs)
// ? list(x)
// : x <= head(xs)
// ? pair(x, xs)
// : pair(head(xs), insert(x, tail(xs)));
// }
// function insertion_sort(xs) {
// return is_null(xs)
// ? xs
// : insert(head(xs), insertion_sort(tail(xs)));
// }
// function list_to_array(L) {
// const A = [];
// let i = 0;
// for (let p = L; !is_null(p); p = tail(p)) {
// A[i] = head(p);
// i = i + 1;
// }
// return A;
// }
// function make_balanced_BST(L) {
// const sorted_L = insertion_sort(L);
// const A = list_to_array(sorted_L);
// function make_tree(low, high) {
// if (low > high) {
// return null;
// } else {
// const mid = math_ceil((low + high) / 2);
// return list(A[mid],
// make_tree(low, mid - 1),
// make_tree(mid + 1, high));
// }
// }
// return make_tree(0, array_length(A) - 1);
// }
// //=====================================================================
// // TASK 3C
// function insert(x, xs) {
// return is_null(xs)
// ? list(x)
// : x <= head(xs)
// ? pair(x, xs)
// : pair(head(xs), insert(x, tail(xs)));
// }
// function insertion_sort(xs) {
// return is_null(xs)
// ? xs
// : insert(head(xs), insertion_sort(tail(xs)));
// }
// // function flatten_bin_tree(T) {
// // return is_null(T)
// // ? null
// // : append(flatten_bin_tree(list_ref(T, 1)),
// // pair(head(T), flatten_bin_tree(list_ref(T, 2))));
// // }
// function bin_tree_to_BST(T) {
// const L = flatten_bin_tree(T);
// let SL = insertion_sort(L);
// function copy_btree(btree) {
// if (is_null(btree)) {
// return null;
// } else {
// const left = copy_btree(list_ref(btree, 1));
// const entry = head(SL);
// SL = tail(SL);
// const right = copy_btree(list_ref(btree, 2));
// return list(entry, left, right);
// }
// }
// return copy_btree(T);
// }
// //=====================================================================