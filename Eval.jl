using IntervalArithmetic
using DataStructures

#Function to print file name for each file read for visual clarity
function printTitle(title)
    n = length(title)+1
    println(repeat(string("-"), n))
    println("$(title):")
    println(repeat(string("-"), n))

end

# Function to generate the polynomial string for visualization
function polynomial_to_string(coeffs)
    terms = []
    n = length(coeffs) - 1

    for(i , coeff) in enumerate(coeffs)
        if coeff != 0
            power = n - i + 1
            term = if power == 0
                "$(abs(coeff))"
            elseif power == 1
                abs(coeff) == 1 ? "x" : "$(abs(coeff))x"
            else
                abs(coeff) == 1 ? "x^$power" : "$(abs(coeff))x^$power"
            end

            # Determine the sign to prepend
            if coeff > 0 && !isempty(terms)
                push!(terms, "+ $term")
            elseif coeff < 0
                push!(terms, "- $term")
            else
                push!(terms, term)
            end
        end
    end

    return join(terms, " ")
end

# Define the polynomial p(x) and its derivative p'(x)
# Uses Horner's Rule
function eval_poly(coeffs, x)
    result = coeffs[1]
    deriv_res = 0
    n = length(coeffs)
    for i in 2:n
        deriv_res = deriv_res * x + result
        result = result * x + coeffs[i]
    end
    return result, deriv_res
end



# Define the predicates
function C(p_interval)
    return !(in_interval(0.0, p_interval)) 
end

# Function to check if endpoints change signs
function endpoints_change_signs(f, interval)
    a = inf(interval)
    b = sup(interval)
    f_a, f_b = f(a), f(b)
    return (f_a * f_b) < 0
end

# Function to check if endpoints contains roots
function endpoints_zero(f, interval)
    a = inf(interval)
    b = sup(interval)
    f_a, f_b = f(a), f(b)
    return (f_a * f_b) == 0
end


# Implement the EVAL algorithm
function eval_root(coeffs, interval, max_iter=100)
    #Store the root intervals found
    roots = Deque{Interval{Float64}}()
    #Stores intermediate bisected intervals 
    queue = Deque{Interval{Float64}}()
    #Stores edge cases where the root lies on the boundary of an interval
    rep_queue = Deque{Interval{Float64}}()
    #Determines if a root is found 
    rootFound = false
    #Early return for faster function execution
    rootCount = 0

    #initial interval
    push!(queue, interval)
    
    #Max iteration set to 100(can change), after which the function terminates 
    for i in 1:max_iter
        
        if isempty(queue)
            break  # Exit early if queue is empty
        end
        
        #Main loop
        while !isempty(queue)
            #Popping intervals in the order of left to right
            I = popfirst!(queue)

            #Evaluate the polynomial on the interval
            f_I, f_p_I = eval_poly(coeffs, I)
            
            #Predicate Checking
            if !C(f_I)
                if !C(f_p_I) 
                    #Bisection Step
                    I_left, I_right = bisect(I, 0.5)
                    push!(queue, I_left)
                    push!(queue, I_right)
                elseif endpoints_change_signs(x -> eval_poly(coeffs, x)[1], I)
                    #Normal root found procedure
                    push!(roots, I)
                    rootFound = true
                    rootCount += 1 
                    if rootCount >= length(coeffs) - 1
                        return roots # Exit early if queue is empty
                    end
                elseif endpoints_zero(x -> eval_poly(coeffs, x)[1], I)
                    #Repeated root found on boundaries
                    push!(rep_queue, I)
                    if(length(rep_queue) == 2)
                        combinedRootInterval = hull(popfirst!(rep_queue), I)
                        push!(roots, combinedRootInterval)
                        rootFound = true
                        rootCount += 1 
                        if rootCount >= length(coeffs) - 1
                            return roots # Exit early if queue is empty
                        end
                    end 
                end 
            end
        end
    end

    #If after max iterations and no root is found, print failure message
    if !rootFound
        println("Maximum number of iterations reached without convergence")
    end

    #Return result if polynomial contains complex and real roots
    return roots
end


#Iterate each file and then each line of the file and apply EVAL on each polynomial and interval given
for i in 1:length(ARGS)
    file = open(ARGS[i])
    printTitle(ARGS[i])
    for line in eachline(file)
        #Skipping comment and empty lines
        if startswith(line, "#") || isempty(line)
            continue
        end

        #Input processing
        description = ""
        parts = split(line, r",\s*")
        if !all(isdigit, parts[end])
            description = parts[end]
            parts = parts[1:end-1]
        end

        coeffs = parse.(Int, parts[1:end-2])
        initial_interval = interval(parse(Float64, parts[end-1]), parse(Float64, parts[end]))

        #Result printing
        if description != ""
            description = replace(description, "\\n" => "\n")
            println(description)
            println()
        end
        
        println(polynomial_to_string(coeffs))

        roots = eval_root(coeffs,initial_interval)
        for root in roots
            println(replace("Root interval found: $root", r"_\w*$" => ""))
        end

        println()
    end
    close(file)
end



