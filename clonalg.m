

search_space = [-5, 5];
max_gens = 100;
pop_size = 100;
clone_factor = 0.1;
num_rand = 2;
bits_per_param = 16;
mutate_factor = -2.5;

best = search(search_space, max_gens, pop_size, clone_factor, num_rand, bits_per_param);

function best = search(search_space, max_gens, pop_size, clone_factor, num_rand, bits_per_param)
    
    pop = cell(pop_size, 1);
    for i = 1 : pop_size
        pop{i} = containers.Map('bitstring', random_bitstring(length(search_space)*bits_per_param), 'UniformValues', false);
    end
   
    evaluate(pop, search_space, bits_per_param);
    sorted_pop = sort_by_cost(pop);
    best = sorted_pop{1};
    
    fprintf("Generation: %d f(%d)\n", 1, 'cost');
    
     for gen = 2 : max_gens
         clones = clone_and_hypermutate(pop, clone_factor);
         evaluate(clones, search_space, bits_per_param);
         
         for x = 1 : length(clones)
            pop{ length(pop) + 1 } = clones{x};
         end
         pop = sort_by_cost(pop);
         pop = pop(1:pop_size);
         pop = random_insertion(search_space, pop, num_rand, bits_per_param);
         best = pop{1};
         
         fprintf("Generation: %d %s\n", gen, best('cost'));
     end
     
        
    best;
end

function bitstring = random_bitstring(num_bits)
    bitstring = "";
    for i = 1 : num_bits
        bit = "";
       if (rand < 0.5) 
           bit = "1";
       else
           bit = "0";
       end
       bitstring = bitstring + bit;
    end
    
    
    bitstring;
end

function evaluate(pop, search_space, bits_per_param)

    for i = 1 : length(pop)
        p = pop{i};
        p('vector') = decode(p('bitstring'), search_space, bits_per_param);        
        p('cost') = objective_function(p('vector'));     
    end
end

function vector = decode(bitstring, search_space, bits_per_param)
    vector = [];
    
    for x = 1: length(search_space)
        
        
        sum = 0;
        
        if x == 1
            off = 1;
            off_end = bits_per_param;
        else
            off = (x - 1) * bits_per_param + 1;
            off_end =   x * bits_per_param ; 
        end
        
        param = reverse(extractBetween(bitstring, off, off_end));
        params_split = split(param, "");
        

        
        for y = 1 : length(params_split)
            add_num = 0;

            if (params_split(y) == "1")
                add_num = 1;
            end
            sum = sum + (add_num * 2.0^(y - 1));
        end
        

        max_bound = max(search_space);
        min_bound = min(search_space);        
       
        new_val = min_bound + ((max_bound - min_bound)/((2.0^bits_per_param) - 1)) * sum;
        
        vector(length(vector)+1) = new_val;
        
    end
    
    vector;
end


function measure = objective_function(vector)
    measure = 0;
    for v = vector
       measure = measure + (v^2); 
    end    
     measure;
end

function sorted = sort_by_cost(pop)
    
    
    sortable = zeros(length(pop));
    switch_prop = cell(length(pop), 1);
   
    for i = 1: length(pop)
        sortable(i) = pop{i}('cost');
    end
    
    [sorted_vals, original_idx] = sort(sortable);
    
    
    for i = 1: length(pop)
         switch_prop{i} = pop{original_idx(i)};
    end
    
    sorted = switch_prop;
end

function child = point_mutation(bitstring, mutation_rate)
    child = "";
        
        bit_str = split(bitstring, "");
        
        
        for i = 1 : length(bit_str)
            bit = bit_str(i);
            random_number = rand;
            if random_number < mutation_rate
                if bit == "1"
                    bit = 0;
                else
                    bit = 1;
                end
            end
            child = child + bit;
        end
 
    child;
end

function clones = clone_and_hypermutate(pop, clone_factor) 
    
    num_clones = floor(length(pop) * clone_factor);
    clones = cell(num_clones, 1);
    calculate_affinity(pop);
    
    for i = 1 : length(pop)
              
        mutation_rate = calculate_mutation_rate(pop{i}, -2.5);
        
        for x = 1: num_clones
            mutated_dna = point_mutation(pop{x}('bitstring'), mutation_rate);
            clones{x} = containers.Map('bitstring', mutated_dna, 'UniformValues', false);
        end
       
    end
    clones;
end

function calculate_affinity(pop)
    pop = sort_by_cost(pop);
    range = pop{length(pop)}('cost') - pop{1}('cost');
    
    if range == 0
        for i = 1 : length(pop)
           pop{i}('affinity') = 1;
        end
    else
        
        for i = 1 : length(pop)
            pop{i}('affinity') = 1-(pop{i}('cost')/range);
        end
    end
end
 
function rate = calculate_mutation_rate(antibody, mutate_factor)
    
    rate = exp(mutate_factor * antibody('affinity'));
end


function population = random_insertion(search_space, pop, num_rand, bits_per_param)
    population = ''; 
    if num_rand == 0
        population = pop;
    else
        random_entities = cell(num_rand, 1);
        for x = 1 : length(random_entities)
            random_entities{x} = containers.Map('bitstring', random_bitstring(length(search_space)*bits_per_param), 'UniformValues', false); 
        end
        
        evaluate(random_entities, search_space, bits_per_param)
        for x = 1 : length(random_entities)
            pop{ length(pop) + 1 } = random_entities{x};
        end
        
        population = sort_by_cost(pop);
    end
    population;
    
    
end

