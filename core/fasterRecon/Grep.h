#pragma once

#ifndef _GREP_H_
#define _GREP_H_

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
//
//
//  -------------------------------------------------------------
//  ----------  MANIFEST CONSTANTS and MACROS BEGIN HERE --------
//  -------------------------------------------------------------
//  NOTE:  These are used internally and awareness of them is not
//         required in the programmer's interface.
//
//  Flags for anding and oring
#define BIT_CAN_MATCH_NOTHING          0x00
#define BIT_MORE_THAN_NULL             0x01
#define BIT_OPENS_ON_STAR_OR_PLUS      0x04
#define BIT_STAR_OR_PLUS_OK            0x02
//  Character values
#define CHAR_CLOSE_BRACKET    ']'
#define CHAR_CLOSE_PAREN      ')'
#define CHAR_DOLLAR           '$'
#define CHAR_DOUBLE_BSLASH    '\\'
#define CHAR_HOCH             '^'
#define CHAR_MINUS            '-'
#define CHAR_OPEN_BRACKET     '['
#define CHAR_OPEN_PAREN       '('
#define CHAR_OR               '|'
#define CHAR_PERIOD           '.'
#define CHAR_PLUS             '+'
#define CHAR_QUESTION         '?'
#define CHAR_STAR             '*'
#define CHAR_ZERO             '\0'
//  Pleasantly obscure macros
#define MAK_AS_UNSIGNED(p)          ((int)*(unsigned char *)(p))
#define MAK_CHAR_2_OPCODE(p)        (*(p))
#define MAK_CONTROLLING_OPERAND(p)  ((p)+3)
#define MAK_GET_NEXT(p)             (((*((p)+1)&0xFF)<<8)+(*((p)+2)&0xFF))
#define MAK_IS_COMPOUND(c)          ((c)==CHAR_STAR||(c)==CHAR_PLUS||(c)==CHAR_QUESTION)
//  Limitations
#define MAX_STORE             1024
#define MAX_SUB_EXPRESSIONS     10
// Opcodes used internally
#define OP_ANY                3
#define OP_ANYBUT             5
#define OP_ANYOF              4
#define OP_BACK               7
#define OP_BOL                1
#define OP_BRANCH             6
#define OP_CLOSE             30
#define OP_END                0
#define OP_EOL                2
#define OP_EXACT              8
#define OP_NMEMPTY            9
#define OP_OPEN              20
#define OP_PLUS              11
#define OP_STAR              10
#define SIGNATURE_VALUE    0x9C
//
//
//--------------------------------------------------------------
//-------------  Error codes that can occur  -------------------
//--------------------------------------------------------------
//  NOTE:  These are the only string literals referenced within
//         the class.  Modify to suit your language preference
//         or desired level of obscurity.
//
extern char *CGrep_errors[];

//--------------------------------------------------------------
//-------------  Class declaration for CGrep  ------------------
//--------------------------------------------------------------
//  NOTE:  As a transferability (not portability, but transferability)
//         consideration, all methods are kept inline.
//
//
class CGrep {
public:
    // query this element to learn if your regular expression was compiled ok.
    int an_error_has_occured;
    // if the CGrep instance did not like something,
    // you can read all about it here.
    char error_string[  255 ];
public:
    // The single constructor requires a UNIX regular
    // expression.  It is ok to pass a constant.
	CGrep( ) { }
	int match(char *exp, char *string ) {
		set_reg_exp(exp);
		return is_reg_exp_present(string);
	}
	void set_reg_exp(char *exp) {
        Special_chars = "^$.[()|?+*\\";
        an_error_has_occured = 0;
        strcpy( error_string, "" );
        if( exp == NULL ) {
            record_error( 0 );
            return;
        }
        Input_scan_pointer = exp;
        parentheses_count = 1;
        symbol_string_size = 0L;
        code_emit_pointer = &reg_exp_shadow;
        emit_next_byte( SIGNATURE_VALUE );
        int flags;
        if( parse_parens( 0, &flags ) == NULL ) {
            record_error( 1 );
            return;
        }
        Input_scan_pointer = exp;
        parentheses_count = 1;
        code_emit_pointer = symbol_string;
        emit_next_byte( SIGNATURE_VALUE );
        if( parse_parens( 0, &flags ) == NULL ) {
            record_error( 1 );
            return;
        }
        first_char_of_regexp = CHAR_ZERO;
        is_anchored = 0;
        required_substring = NULL;
        length_required_substring = 0;
        char *scan = symbol_string + 1;
        if( MAK_CHAR_2_OPCODE( get_next_pointer( scan ) ) == OP_END ) {
            scan = MAK_CONTROLLING_OPERAND( scan );
            if( MAK_CHAR_2_OPCODE( scan ) == OP_EXACT )
                first_char_of_regexp = *MAK_CONTROLLING_OPERAND( scan );
            else if( MAK_CHAR_2_OPCODE( scan ) == OP_BOL )
                is_anchored++;
            if( flags & BIT_OPENS_ON_STAR_OR_PLUS ) {
                char *longest = NULL;
                int len = 0;
                for( ; scan != NULL; scan = get_next_pointer( scan ) )
                     if( MAK_CHAR_2_OPCODE( scan ) == OP_EXACT && strlen( MAK_CONTROLLING_OPERAND( scan ) ) >= len ) {
                         longest = MAK_CONTROLLING_OPERAND( scan );
                         len = strlen( MAK_CONTROLLING_OPERAND( scan ) );
                     }
                required_substring = longest;
                length_required_substring = len;
            }
        }
    }
    //  This method returns 1 if the regular expression
    //  was matched.  Otherwise it returns 0.
    int is_reg_exp_present( char *string ) {
        char *s;
        if( string == NULL ) {
            record_error( 0 );
            return( 0 );
        }
        if( MAK_AS_UNSIGNED( symbol_string ) != SIGNATURE_VALUE ) {
            record_error( 2 );
            return( 0 );
        }
        if( required_substring != NULL ) {
            s = string;
            while( ( s = strchr( s, required_substring[ 0 ] ) ) != NULL ) {
                   if( strncmp( s, required_substring, length_required_substring ) == 0 )
                       break;
                   s++;
            }
            if( s == NULL )
                return( 0 );
        }
        beginning_of_input = string;
        if( is_anchored )
            return( evaluate_next( string ) );
        s = string;
        if( first_char_of_regexp != CHAR_ZERO )
            while( ( s = strchr( s, first_char_of_regexp ) ) != NULL ) {
                   if( evaluate_next( s ) )
                       return( 1 );
                   s++;
            }
        else do {
            if( evaluate_next( s ) )
                return( 1 );
        } while( *s++ != CHAR_ZERO );
        return( 0 );
    }
private:
    //----------------------------------------------------------
    //------------  THE PRIVATE INTERFACE BEGINS HERE  ---------
    //----------------------------------------------------------
    char *reg_exp_match_begins_here[ MAX_SUB_EXPRESSIONS ];
    char *reg_exp_match_ends_here[ MAX_SUB_EXPRESSIONS ];
    char first_char_of_regexp;
    char is_anchored;
    char *required_substring;
    int length_required_substring;
    char symbol_string[ MAX_STORE ];
    char *Input_scan_pointer;
    int parentheses_count;
    char reg_exp_shadow;
    char *code_emit_pointer;
    long symbol_string_size;
    char *input_string_pointer;
    char *beginning_of_input;
    char **pointer_to_rembh;
    char **pointer_to_remeh;
    char *Special_chars;
    //----------------------------------------------------------
    //----------------------------------------------------------
    void record_error( int n ) {
        strcpy( error_string, CGrep_errors[ n ] );
        ++an_error_has_occured;
    }
    //----------------------------------------------------------
    char *parse_parens( int paren, int *flagp ) {
        *flagp = BIT_MORE_THAN_NULL;
        char *ret = NULL;
        int parno = 0;
        if( paren ) {
            if( parentheses_count >= MAX_SUB_EXPRESSIONS ) {
                record_error( 3 );
                return( NULL );
            }
            parno = parentheses_count;
            parentheses_count++;
            ret = emit_node( OP_OPEN + parno );
        } else
            ret = NULL;
        int flags;
        char *br = one_side_of_or_operator( &flags );
        if( br == NULL )
            return( NULL );
        if( ret != NULL )
            go_to_end_of_chain( ret, br );
        else
            ret = br;
        if( !( flags & BIT_MORE_THAN_NULL ) )
            *flagp &= ~BIT_MORE_THAN_NULL;
        *flagp |= flags & BIT_OPENS_ON_STAR_OR_PLUS;
        while( *Input_scan_pointer == CHAR_OR ) {
               Input_scan_pointer++;
               br = one_side_of_or_operator( &flags );
               if( br == NULL )
                   return( NULL );
               go_to_end_of_chain( ret, br );
               if( !( flags&BIT_MORE_THAN_NULL ) )
                   *flagp &= ~BIT_MORE_THAN_NULL;
                   *flagp |= flags & BIT_OPENS_ON_STAR_OR_PLUS;
        }
        char *ender = emit_node( ( paren ) ? OP_CLOSE + parno : OP_END );
        go_to_end_of_chain( ret, ender );
        for( br = ret; br != NULL; br = get_next_pointer( br ) )
             go_to_end_of_chain_on_operand( br, ender );
        if( paren && *Input_scan_pointer++ != CHAR_CLOSE_PAREN ) {
            record_error( 3 );
            return( NULL );
        } else if( !paren && *Input_scan_pointer != CHAR_ZERO ) {
            if( *Input_scan_pointer == CHAR_CLOSE_PAREN ) {
                record_error( 3 );
                return( NULL );
            } else
                record_error( 4 );
                return( NULL );
        }
        return( ret );
    }
    //----------------------------------------------------------
    char *one_side_of_or_operator( int *flagp ) {
        int flags = 0;
        char *ret = NULL;
        *flagp = BIT_CAN_MATCH_NOTHING;
        ret = emit_node( OP_BRANCH );
        char *chain = NULL;
        char *latest = NULL;
        while( *Input_scan_pointer != CHAR_ZERO && *Input_scan_pointer != CHAR_OR && *Input_scan_pointer != CHAR_CLOSE_PAREN ) {
               latest = trailing_wild( &flags );
               if( latest == NULL )
                   return( NULL );
               *flagp |= flags & BIT_MORE_THAN_NULL;
               if( chain == NULL ) /* First piece. */
                   *flagp |= flags & BIT_OPENS_ON_STAR_OR_PLUS;
               else
                   go_to_end_of_chain( chain, latest );
               chain = latest;
        }
        if( chain == NULL ) /* Loop ran zero times. */
            (void) emit_node( OP_NMEMPTY );
        return( ret );
    }
    //----------------------------------------------------------
    char *trailing_wild( int *flagp ) {
        int flags = 0;
        char *ret = compress_ordinary_characters( &flags );
        char *next = NULL;
        if( ret == NULL )
            return( NULL );
        char op = *Input_scan_pointer;
        if( !MAK_IS_COMPOUND( op ) ) {
            *flagp = flags;
            return( ret );
        }
        if( !( flags & BIT_MORE_THAN_NULL ) && op != CHAR_QUESTION ) {
            record_error( 5 );
            return( NULL );
        }
        *flagp =( op != CHAR_PLUS ) ?( BIT_CAN_MATCH_NOTHING | BIT_OPENS_ON_STAR_OR_PLUS ) :
                                     ( BIT_CAN_MATCH_NOTHING | BIT_MORE_THAN_NULL );
        if( op == CHAR_STAR &&( flags & BIT_STAR_OR_PLUS_OK ) )
            insert_operator( OP_STAR, ret );
        else if( op == CHAR_STAR ) {
            insert_operator( OP_BRANCH, ret );
            go_to_end_of_chain_on_operand( ret, emit_node( OP_BACK ) );
            go_to_end_of_chain_on_operand( ret, ret );
            go_to_end_of_chain( ret, emit_node( OP_BRANCH ) );
            go_to_end_of_chain( ret, emit_node( OP_NMEMPTY ) );
        } else if( op == CHAR_PLUS &&( flags&BIT_STAR_OR_PLUS_OK ) )
            insert_operator( OP_PLUS, ret );
        else if( op == CHAR_PLUS ) {
            next = emit_node( OP_BRANCH );
            go_to_end_of_chain( ret, next );
            go_to_end_of_chain( emit_node( OP_BACK ), ret );
            go_to_end_of_chain( next, emit_node( OP_BRANCH ) );
            go_to_end_of_chain( ret, emit_node( OP_NMEMPTY ) );
        } else if( op == CHAR_QUESTION ) {
            insert_operator( OP_BRANCH, ret );
            go_to_end_of_chain( ret, emit_node( OP_BRANCH ) );
            next = emit_node( OP_NMEMPTY );
            go_to_end_of_chain( ret, next );
            go_to_end_of_chain_on_operand( ret, next );
        }
        Input_scan_pointer++;
        if( MAK_IS_COMPOUND( *Input_scan_pointer ) ) {
            record_error( 6 );
            return( NULL );
        }
        return( ret );
    }
    //----------------------------------------------------------
    char *compress_ordinary_characters( int *flagp ) {
       char *ret = NULL;
       int flags = 0;
       *flagp = BIT_CAN_MATCH_NOTHING;
       switch( *Input_scan_pointer++ ) {
               case CHAR_HOCH:
                    ret = emit_node( OP_BOL );
                    break;
               case CHAR_DOLLAR:
                    ret = emit_node( OP_EOL );
                    break;
               case CHAR_PERIOD:
                    ret = emit_node( OP_ANY );
                    *flagp |= BIT_MORE_THAN_NULL | BIT_STAR_OR_PLUS_OK;
                    break;
               case CHAR_OPEN_BRACKET: {
                    int gjv_reg_class;
                    int classend;
                    if( *Input_scan_pointer == CHAR_HOCH ) {
                        ret = emit_node( OP_ANYBUT );
                        Input_scan_pointer++;
                    } else
                        ret = emit_node( OP_ANYOF );
                    if( *Input_scan_pointer == CHAR_CLOSE_BRACKET || *Input_scan_pointer == CHAR_MINUS )
                        emit_next_byte( *Input_scan_pointer++ );
                    while( *Input_scan_pointer != CHAR_ZERO && *Input_scan_pointer != CHAR_CLOSE_BRACKET ) {
                           if( *Input_scan_pointer == CHAR_MINUS ) {
                               Input_scan_pointer++;
                               if( *Input_scan_pointer == CHAR_CLOSE_BRACKET || *Input_scan_pointer == CHAR_ZERO )
                                   emit_next_byte( CHAR_MINUS );
                               else {
                                   gjv_reg_class = MAK_AS_UNSIGNED( Input_scan_pointer - 2 ) + 1;
                                   classend = MAK_AS_UNSIGNED( Input_scan_pointer );
                                   if( gjv_reg_class > classend + 1 ) {
                                       record_error( 7 );
                                       return( NULL );
                                   }
                                   for( ; gjv_reg_class <= classend; gjv_reg_class++ )
                                        emit_next_byte( gjv_reg_class );
                                   Input_scan_pointer++;
                               }
                           } else
                               emit_next_byte( *Input_scan_pointer++ );
                    }
                    emit_next_byte( CHAR_ZERO );
                    if( *Input_scan_pointer != CHAR_CLOSE_BRACKET ) {
                        record_error( 3 );
                        return( NULL );
                    }
                    Input_scan_pointer++;
                    *flagp |= BIT_MORE_THAN_NULL | BIT_STAR_OR_PLUS_OK; }
                    break;
               case CHAR_OPEN_PAREN:
                    ret = parse_parens( 1, &flags );
                    if( ret == NULL )
                        return( NULL );
                    *flagp |= flags & ( BIT_MORE_THAN_NULL | BIT_OPENS_ON_STAR_OR_PLUS );
                    break;
               case CHAR_ZERO:
               case CHAR_OR:
               case CHAR_CLOSE_PAREN:
                    record_error( 2 );
                    return( NULL );
               case CHAR_QUESTION:
               case CHAR_PLUS:
               case CHAR_STAR:
                    record_error( 8 );
                    return( NULL );
               case CHAR_DOUBLE_BSLASH:
                    if( *Input_scan_pointer == CHAR_ZERO ) {
                        record_error( 9 );
                        return( NULL );
                    }
                    ret = emit_node( OP_EXACT );
                    emit_next_byte( *Input_scan_pointer++ );
                    emit_next_byte( CHAR_ZERO );
                    *flagp |= BIT_MORE_THAN_NULL | BIT_STAR_OR_PLUS_OK;
                    break;
               default: {
                    int len;
                    char ender;
                    Input_scan_pointer--;
                    len = strcspn( Input_scan_pointer, Special_chars );
                    if( len <= 0 ) {
                        record_error( 2 );
                        return( NULL );
                    }
                    ender = *( Input_scan_pointer+len );
                    if( len > 1 && MAK_IS_COMPOUND( ender ) )
                        len--;
                    *flagp |= BIT_MORE_THAN_NULL;
                    if( len == 1 )
                        *flagp |= BIT_STAR_OR_PLUS_OK;
                    ret = emit_node( OP_EXACT );
                    while( len > 0 ) {
                           emit_next_byte( *Input_scan_pointer++ );
                           len--;
                    }
                    emit_next_byte( CHAR_ZERO ); }
                    break;
       }
       return( ret );
    }
    //----------------------------------------------------------
    char *emit_node( char op ) {
        char *ret = code_emit_pointer;
        if( ret == &reg_exp_shadow ) {
            symbol_string_size += 3;
            return( ret );
        }
        char *ptr = ret;
        *ptr++ = op;
        *ptr++ = CHAR_ZERO;
        *ptr++ = CHAR_ZERO;
        code_emit_pointer = ptr;
        return( ret );
    }
    //----------------------------------------------------------
    void emit_next_byte( char b ) {
        if( code_emit_pointer != &reg_exp_shadow )
            *code_emit_pointer++ = b;
        else
            symbol_string_size++;
    }
    //----------------------------------------------------------
    void insert_operator( char op, char *opnd ) {
        if( code_emit_pointer == &reg_exp_shadow ) {
            symbol_string_size += 3;
            return;
        }
        char *src = code_emit_pointer;
        code_emit_pointer += 3;
        char *dst = code_emit_pointer;
        while( src > opnd )
               *--dst = *--src;
        char *place = opnd;
        *place++ = op;
        *place++ = CHAR_ZERO;
        *place++ = CHAR_ZERO;
    }
    //----------------------------------------------------------
    void go_to_end_of_chain( char *p, char *val ) {
        int offset = 0;
        if( p == &reg_exp_shadow )
            return;
        char *scan = p;
        for( ;; ) {
             char *temp = get_next_pointer( scan );
             if( temp == NULL )
                 break;
             scan = temp;
        }
        if( MAK_CHAR_2_OPCODE( scan ) == OP_BACK )
            offset = scan - val;
        else
            offset = val - scan;
        *( scan + 1 ) = ( offset >> 8 ) & 0xFF;
        *( scan + 2 ) = offset & 0xFF;
    }
    //----------------------------------------------------------
    void go_to_end_of_chain_on_operand( char *p, char *val ) {
        if( p == NULL || p == &reg_exp_shadow || MAK_CHAR_2_OPCODE( p ) != OP_BRANCH )
            return;
        go_to_end_of_chain( MAK_CONTROLLING_OPERAND( p ), val );
    }
    //----------------------------------------------------------
    int lookup_engine( char *prog ) {
        char *scan = prog;
        char *next = NULL;
        while( scan != NULL ) {
               next = get_next_pointer( scan );
               switch( MAK_CHAR_2_OPCODE( scan ) ) {
                       case OP_BOL:
                            if( input_string_pointer != beginning_of_input )
                                return( 0 );
                            break;
                       case OP_EOL:
                            if( *input_string_pointer != CHAR_ZERO )
                                return( 0 );
                            break;
                       case OP_ANY:
                            if( *input_string_pointer == CHAR_ZERO )
                                return( 0 );
                            input_string_pointer++;
                            break;
                       case OP_EXACT: {
                            int len;
                            char *opnd;
                            opnd = MAK_CONTROLLING_OPERAND( scan );
                            if( *opnd != *input_string_pointer )
                                return( 0 );
                            len = strlen( opnd );
                            if( len > 1 && strncmp( opnd, input_string_pointer, len ) != 0 )
                                return( 0 );
                            input_string_pointer += len;}
                            break;
                       case OP_ANYOF:
                            if( *input_string_pointer == CHAR_ZERO || strchr( MAK_CONTROLLING_OPERAND( scan ), *input_string_pointer ) == NULL )
                                return( 0 );
                            input_string_pointer++;
                            break;
                       case OP_ANYBUT:
                            if( *input_string_pointer == CHAR_ZERO || strchr( MAK_CONTROLLING_OPERAND( scan ), *input_string_pointer ) != NULL )
                                return( 0 );
                            input_string_pointer++;
                            break;
                       case OP_NMEMPTY:
                            break;
                       case OP_BACK:
                            break;
                       case OP_OPEN + 1:
                       case OP_OPEN + 2:
                       case OP_OPEN + 3:
                       case OP_OPEN + 4:
                       case OP_OPEN + 5:
                       case OP_OPEN + 6:
                       case OP_OPEN + 7:
                       case OP_OPEN + 8:
                       case OP_OPEN + 9: {
                            int no;
                            char *save;
                            no = MAK_CHAR_2_OPCODE( scan ) - OP_OPEN;
                            save = input_string_pointer;
                            if( lookup_engine( next ) ) {
                                if( pointer_to_rembh[ no ] == NULL )
                                    pointer_to_rembh[ no ] = save;
                                return( 1 );
                            } else
                                return( 0 );}
                            break;
                       case OP_CLOSE + 1:
                       case OP_CLOSE + 2:
                       case OP_CLOSE + 3:
                       case OP_CLOSE + 4:
                       case OP_CLOSE + 5:
                       case OP_CLOSE + 6:
                       case OP_CLOSE + 7:
                       case OP_CLOSE + 8:
                       case OP_CLOSE + 9: {
                            int no;
                            char *save;
                            no = MAK_CHAR_2_OPCODE( scan ) - OP_CLOSE;
                            save = input_string_pointer;
                            if( lookup_engine( next ) ) {
                                if( pointer_to_remeh[ no ] == NULL )
                                    pointer_to_remeh[ no ] = save;
                                return( 1 );
                            } else
                                return( 0 ); }
                            break;
                       case OP_BRANCH: {
                            char *save;
                            if( MAK_CHAR_2_OPCODE( next ) != OP_BRANCH )
                                next = MAK_CONTROLLING_OPERAND( scan );
                            else {
                                do {
                                   save = input_string_pointer;
                                   if( lookup_engine( MAK_CONTROLLING_OPERAND( scan ) ) )
                                       return( 1 );
                                   input_string_pointer = save;
                                   scan = get_next_pointer( scan );
                                } while( scan != NULL && MAK_CHAR_2_OPCODE( scan ) == OP_BRANCH );
                                return( 0 );
                            }}
                            break;
                       case OP_STAR:
                       case OP_PLUS: {
                            char nextch;
                            int no;
                            char *save;
                            int min;
                            nextch = CHAR_ZERO;
                            if( MAK_CHAR_2_OPCODE( next ) == OP_EXACT )
                                nextch = *MAK_CONTROLLING_OPERAND( next );
                            min = ( MAK_CHAR_2_OPCODE( scan ) == OP_STAR ) ? 0 : 1;
                            save = input_string_pointer;
                            no = wildcard_lookup( MAK_CONTROLLING_OPERAND( scan ) );
                            while( no >= min ) {
                                   if( nextch == CHAR_ZERO || *input_string_pointer == nextch )
                                       if( lookup_engine( next ) )
                                           return( 1 );
                                   no--;
                                   input_string_pointer = save + no;
                            }
                            return( 0 );}
                            break;
                       case OP_END:
                            return( 1 );
                       default:
                            record_error( 2 );
                            return( 0 );
               }
               scan = next;
        }
        record_error( 2 );
        return( 0 );
    }
    //----------------------------------------------------------
    int wildcard_lookup( char *p ) {
        int count = 0;
        char *scan = input_string_pointer;
        char *opnd = MAK_CONTROLLING_OPERAND( p );
        switch( MAK_CHAR_2_OPCODE( p ) ) {
                case OP_ANY:
                     count = strlen( scan );
                     scan += count;
                     break;
                case OP_EXACT:
                     while( *opnd == *scan ) {
                            count++;
                            scan++;
                     }
                     break;
                case OP_ANYOF:
                     while( *scan != CHAR_ZERO && strchr( opnd, *scan ) != NULL ) {
                            count++;
                            scan++;
                     }
                     break;
                case OP_ANYBUT:
                     while( *scan != CHAR_ZERO && strchr( opnd, *scan ) == NULL ) {
                            count++;
                            scan++;
                     }
                     break;
                default:
                     record_error( 2 );
                     count = 0;
                     return( NULL );
        }
        input_string_pointer = scan;
        return( count );
    }
    //----------------------------------------------------------
    char *get_next_pointer( char *p ) {
        if( p == &reg_exp_shadow )
            return( NULL );
        int offset = MAK_GET_NEXT( p );
        if( offset == 0 )
            return( NULL );
        if( MAK_CHAR_2_OPCODE( p ) == OP_BACK )
            return( p - offset );
        else
            return( p + offset );
    }
    //----------------------------------------------------------
    int evaluate_next( char *string ) {
        input_string_pointer = string;
        pointer_to_rembh = reg_exp_match_begins_here;
        pointer_to_remeh = reg_exp_match_ends_here;
        char **sp = reg_exp_match_begins_here;
        char **ep = reg_exp_match_ends_here;
        for( int i = MAX_SUB_EXPRESSIONS; i > 0; i-- ) {
             *sp++ = NULL;
             *ep++ = NULL;
        }
        if( lookup_engine( symbol_string + 1 ) ) {
            reg_exp_match_begins_here[ 0 ] = string;
            reg_exp_match_ends_here[ 0 ] = input_string_pointer;
            return( 1 );
        } else
           return( 0 );
    }
};

float get_radius(char* resname, char* aname);

#endif