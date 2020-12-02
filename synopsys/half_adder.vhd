LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.numeric_std.all;

ENTITY half_adder IS
PORT
 (		
		A : IN STD_LOGIC;
		B : IN STD_LOGIC;
		S : OUT STD_LOGIC;
		C : OUT STD_LOGIC
 );
END half_adder;

ARCHITECTURE Behavior OF half_adder IS
BEGIN

	S <= A XOR B;
	C <= A AND B;

END Behavior;