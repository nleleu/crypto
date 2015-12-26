require 'prime'
n = 1073

Prime.each do |p|

puts("Test : #{p}")

if (n%p)==0
puts ("Prime1 #{p}")
puts ("Prime2 #{n/p}")
break
end
end
