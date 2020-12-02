LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.numeric_std.all;

ENTITY full_adder IS
PORT
 (		
		A 			: IN STD_LOGIC;
		B 			: IN STD_LOGIC;
		carry_in	: IN STD_LOGIC;
		S 			: OUT STD_LOGIC;
		C 			: OUT STD_LOGIC
 );
END full_adder;

ARCHITECTURE Behavior OF full_adder IS

COMPONENT half_adder IS
PORT
 (		
		A : IN STD_LOGIC;
		B : IN STD_LOGIC;
		S : OUT STD_LOGIC;
		C : OUT STD_LOGIC
 );
END COMPONENT;

SIGNAL ha1_s, ha1_c, ha2_c 	:	STD_LOGIC;


BEGIN


Ha1: half_adder
PORT MAP(A, B, ha1_s, ha1_c);

Ha2: half_adder
PORT MAP(ha1_s, carry_in, S, ha2_c);

C <= ha1_c OR ha2_c;


END Behavior;