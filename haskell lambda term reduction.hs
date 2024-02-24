--typ
data LambdaTerm = Var Char
  | App LambdaTerm LambdaTerm
  | Lam Char LambdaTerm
    deriving (Eq, Show,Read)


--jakies przyklady
exampleExpr1 :: LambdaTerm
exampleExpr1 = Lam 'z' (App (Lam 'x' ( App (Var 'x') (  Var 'x') )) (Var 'y'))


exampleExpr2 :: LambdaTerm
exampleExpr2 = Lam 's' (Lam 'x'  (App (Var 's')  (App (Var 's') ( App (Var 's')  (App (Var 's') (Var 'x'))))))

exampleExpr3 :: LambdaTerm
exampleExpr3 = Lam 'y' (Lam 'z'  (App (Var 'y')  (App (Var 'y') ( App (Var 'y')  (App (Var 'y') (Var 'z'))))))


exampleExpr4 :: LambdaTerm
exampleExpr4 = Lam 's' (Lam 'x'  (App (Var 's')  (App  (Var 's')  (Var 'x') )    ))


exampleExpr5 :: LambdaTerm
exampleExpr5 = Lam 't' (Lam 'y'  (App (Var 't')  (App  (Var 't')  (Var 'y') )    ))

--funkcja do wizaulizacji
visualize :: LambdaTerm -> String
visualize (Var x) = [x]
visualize (App e1 e2) = "(" ++ visualize e1 ++ " " ++ visualize e2 ++ ")"
visualize (Lam x body) = "L" ++ [x] ++ "." ++ visualize body


--funkcja do podmiany w betaredukcji
substitute :: LambdaTerm -> Char -> LambdaTerm -> LambdaTerm
substitute (Var v) x replacement
  | v == x    = replacement
  | otherwise = Var v
substitute (App t1 t2) x replacement = App (substitute t1 x replacement) (substitute t2 x replacement)
substitute (Lam v body) x replacement
  | v == x    = Lam v body
  | otherwise = Lam v (substitute body x replacement)


--zwraca zmienne w wyrazeniu
getVariables :: LambdaTerm -> [Char]
getVariables (Var v)      = [v]
getVariables (App t1 t2)  = getVariables t1 ++ getVariables t2
getVariables (Lam v body) = v : getVariables body

--alfakonwersja
alpha :: Char -> Char -> LambdaTerm->LambdaTerm
alpha o n (Var v)      | v==o = (Var n)
                       | True = (Var v)
alpha o n (App t1 t2)  = (App (alpha o n t1) (alpha o n t2))
alpha o n (Lam v body) | v==o = (Lam n (alpha o n body ))
                       | True = (Lam v (alpha o n body ))

--to jest mi potrzebne dalej
alphaUC = uncurry alpha

--sprawdza czy zmienna o danym znaku jest wolna
isFreeIn :: Char -> LambdaTerm -> Bool
isFreeIn x (Var v)      = x == v
isFreeIn x (App t1 t2)  = isFreeIn x t1 || isFreeIn x t2
isFreeIn x (Lam v body) = x /= v && isFreeIn x body


--to jest mi potrzebne dalej, fukcja zwraca tylko odwrotna wartosc do isfreein
isBoundIn :: Char -> LambdaTerm -> Bool
isBoundIn x y = not (isFreeIn x y)

--usuwa znak ze stringa
removeChar :: String -> Char -> String
removeChar "" c = ""
removeChar (a:b) c | (a==c) = (removeChar b c)
                   | True = a: (removeChar b c)

--usuwa znaki ktore sa w stringu 2 ze stringa 1 
removeString :: String -> String -> String
removeString "" x = ""
removeString x "" = x
removeString a (c:d) = removeString (removeChar a c) d 

--generuje mi nowe zmienne ktore nie sa w uzyciu
getListNew :: String -> Int -> String
getListNew x n = take n (removeString ['a'..'z'] x )


--usuwa duplikaty ze stringu
removeduplicates :: String-> String
removeduplicates "" = ""
removeduplicates (a:b) = [a] ++ (removeduplicates (removeChar b a)) 

--zwraca mi nowe mapowanie jak alfa konwertowac lewy lambda term zeby nie doszlo do kolizji
preventColisionTable :: LambdaTerm -> LambdaTerm -> [(Char,Char)]
preventColisionTable termL termR = zip ( filter ((flip isBoundIn) termL)   (removeduplicates (getVariables termL))  ) (getListNew (  (getVariables termR)++(getVariables termL)) (length (getVariables termL)))

--zamienia lewy lambda term tak by nie bylo kolizji
preventColision :: LambdaTerm -> LambdaTerm -> LambdaTerm
preventColision termL termR = foldl (flip alphaUC) (termL) (preventColisionTable termL termR)


--betaredukcja naiwna nie jest dalej uzyta
betaReduce :: LambdaTerm -> LambdaTerm
betaReduce (App (Lam x body) arg) = substitute body x arg
betaReduce term                   = term

--betaredukcja w normalizacji zamienia zmienne w lewym termie przed aplikacja tak by nie bylo kolizji
betaReduceSafe :: LambdaTerm -> LambdaTerm
betaReduceSafe (App (Lam x body) arg) = substitute (preventColision body arg   ) x arg
betaReduceSafe term                   = term


--funkcja normalizujaca term ktory jest aplikajca
normalizeH :: LambdaTerm -> LambdaTerm
normalizeH term =
  let reduced = betaReduceSafe term
  in if reduced == term then term else normalizeH reduced


--funkcja rekurencyjnie normalizujaca term
normalize :: LambdaTerm -> LambdaTerm
normalize (Var x) = (Var x)
normalize (Lam x body) = (Lam x (normalize body))
normalize (App e1 e2) | (normalizeH (App e1 e2)) /= (App e1 e2) = normalize (normalizeH (App e1 e2))
                      | ((App (normalize e1) e2)) /= (App e1 e2) = normalize ((App (normalize e1) e2))
                      | ((App e1 (normalize e2))) /= (App e1 e2) = normalize ((App e1 (normalize e2)))
                      | True = (App e1 e2)


replacePatterns :: String -> String
replacePatterns [] = []
replacePatterns ('/' : c : '.' : rest) = "Lam '" ++ [c] ++ "' " ++ replacePatterns rest
replacePatterns ('(' : c : ')' : rest) = "(Var '" ++ [c] ++ "')" ++ replacePatterns rest
replacePatterns (x : xs) = x : replacePatterns xs

replaceDoubleParentheses :: String -> String
replaceDoubleParentheses str =
  let updatedStr = replaceOnce str
  in if updatedStr == str
     then str
     else replaceDoubleParentheses updatedStr

replaceOnce :: String -> String
replaceOnce [] = []
replaceOnce ('(': '(' : rest) = "( App (" ++ replaceOnce rest
replaceOnce (x:xs) = x : replaceOnce xs


main :: IO ()
main = do
  putStrLn "Podaj term lambda:"
  input <- getLine  
  let inputC = replaceDoubleParentheses (replacePatterns (input))
  putStrLn inputC
  let myTerm = (read (inputC) :: LambdaTerm)
  putStrLn ( "Przyjety term lambda to       :   " ++ (show  (visualize(myTerm))  ))
  putStrLn ( "Term lambda po normalizacji to:   " ++ (show  (visualize(normalize myTerm ))))






