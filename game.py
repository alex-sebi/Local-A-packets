from pa import *
from EMS import *
from graph import *
from pic1 import *
from draw import *




import pygame
import sys

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 600
BOX_HEIGHT = 30
VERTICAL_SPACING = 15
BOTTOM_MARGIN = 100
LEFT_MARGIN = 100
RIGHT_MARGIN = 200
BOX_OFFSET_Y = 40

DIFF_BOX_WIDTH = 60
DIFF_BOX_HEIGHT = 30
GAP_BETWEEN_BOXES = 10

BUTTON_WIDTH = 140
BUTTON_HEIGHT = 35
BUTTON_MARGIN = 10

# Colors
WHITE = (255, 255, 255)
LIGHT_BLUE = (173, 216, 230)
RED = (200, 0, 0)
GREEN = (0, 200, 0)
BLACK = (0, 0, 0)
GRAY = (220, 220, 220)
BLUE = (0, 100, 255)
DARK_GRAY = (100, 100, 100)
BUTTON_BG = (180, 180, 255)
BUTTON_BG_HOVER = (150, 150, 255)

def convert_to_int_tuples(M):
    return [([int(A), int(B)], int(m)) for ([A, B], m) in M]


def sw(labels,i):
    l=[]
    for k in range(0,i):
        l.append(labels[k])
    l.append(labels[i+1])
    l.append(labels[i])
    for k in range(i+2,len(labels)):
        l.append(labels[k])
    return(l)


def draw_main_view(screen, M, labels, font, grid_start, grid_end, grid_spacing, button_rects):
    screen.fill(WHITE)

    for (label, rect) in button_rects:
        mouse_pos = pygame.mouse.get_pos()
        is_hover = rect.collidepoint(mouse_pos)
        bg = BUTTON_BG_HOVER if is_hover else BUTTON_BG
        pygame.draw.rect(screen, bg, rect)
        pygame.draw.rect(screen, BLACK, rect, 1)
        label_surf = font.render(label, True, BLACK)
        screen.blit(label_surf, label_surf.get_rect(center=rect.center))

    def coord_to_x(x_val):
        return LEFT_MARGIN + (x_val - grid_start) * grid_spacing

    num_items = len(M)
    box_centers = []  
    box_tops = []     

    # Draw axes & labels
    for x in range(grid_start, grid_end + 1):
        x_pos = coord_to_x(x)
        pygame.draw.line(screen, GRAY, (x_pos+grid_spacing/2, 0 + BOX_OFFSET_Y), (x_pos+grid_spacing/2, SCREEN_HEIGHT - BOTTOM_MARGIN), 1)
        # Tick
        pygame.draw.line(screen, BLACK, (x_pos, SCREEN_HEIGHT - BOTTOM_MARGIN),
                         (x_pos, SCREEN_HEIGHT - BOTTOM_MARGIN + 5), 2)
        label = font.render(str(x), True, BLACK)
        screen.blit(label, label.get_rect(center=(x_pos, SCREEN_HEIGHT - BOTTOM_MARGIN + 20)))

    # Draw main boxes
    for idx, ((A, B), m) in enumerate(M):
        x_start = coord_to_x(B)-grid_spacing/2
        x_end = coord_to_x(A)+grid_spacing/2
        full_width = x_end - x_start

        y = (SCREEN_HEIGHT - BOTTOM_MARGIN
             - (BOX_HEIGHT + VERTICAL_SPACING) * (idx)
             - BOX_OFFSET_Y)
        box_tops.append(y)
        center = y + BOX_HEIGHT / 2
        box_centers.append(center)

                # Number label 
        number_label = font.render(str(labels[idx]), True, BLACK)
        label_x = LEFT_MARGIN - 30  
        label_y = y + BOX_HEIGHT / 2
        screen.blit(number_label, number_label.get_rect(center=(label_x, label_y)))

        pygame.draw.rect(screen, LIGHT_BLUE, (x_start, y, full_width, BOX_HEIGHT))
        m_abs = abs(m)
        center_width = m_abs * grid_spacing
        center_color = GREEN if m > 0 else RED
        center_x = x_start + (full_width - center_width) / 2

        pygame.draw.rect(screen, center_color, (center_x, y, center_width, BOX_HEIGHT))
        m_label = font.render(str(m), True, WHITE)
        screen.blit(m_label, m_label.get_rect(center=(center_x + center_width / 2, y + BOX_HEIGHT / 2)))
        pygame.draw.rect(screen, BLACK, (x_start, y, full_width, BOX_HEIGHT), 1)

    # Draw difference & distance clickable boxes
    diff_boxes = []
    distance_boxes = []

    for i in range(len(M) - 1):
        c1 = box_centers[i]
        c2 = box_centers[i + 1]
        mid_center = (c1 + c2) / 2
        diff = M[i + 1][1] - M[i][1]

        # difference box
        diff_color = GREEN if diff > 0 else RED
        diff_x = SCREEN_WIDTH - RIGHT_MARGIN + 30
        diff_y = mid_center - (DIFF_BOX_HEIGHT + VERTICAL_SPACING)/2 
        diff_rect = pygame.Rect(diff_x, diff_y, DIFF_BOX_WIDTH, DIFF_BOX_HEIGHT)
        diff_boxes.append((diff_rect, i))

        pygame.draw.rect(screen, diff_color, diff_rect)
        pygame.draw.rect(screen, BLACK, diff_rect, 1)
        diff_label = font.render(str(diff), True, WHITE)
        screen.blit(diff_label, diff_label.get_rect(center=diff_rect.center))

        # distance box
        (A1, B1), _ = M[i]
        (A2, B2), _ = M[i + 1]
        distance_val = abs(A2 - A1) + abs(B2 - B1)
        dist_x = diff_x + DIFF_BOX_WIDTH + GAP_BETWEEN_BOXES
        dist_y = diff_y
        dist_rect = pygame.Rect(dist_x, dist_y, DIFF_BOX_WIDTH, DIFF_BOX_HEIGHT)

        # clickable only if abs(diff) == distance_val
        if abs(diff) == distance_val:
            distance_boxes.append((dist_rect, i))
            dist_color = BLUE
        else:
            dist_color = DARK_GRAY

        pygame.draw.rect(screen, dist_color, dist_rect)
        pygame.draw.rect(screen, BLACK, dist_rect, 1)
        dist_label = font.render(str(distance_val), True, WHITE)
        screen.blit(dist_label, dist_label.get_rect(center=dist_rect.center))

    return diff_boxes, distance_boxes

def draw_langlands_view(screen, font, button_rects):
    screen.fill(WHITE)
    for (label, rect) in button_rects:
        mouse_pos = pygame.mouse.get_pos()
        is_hover = rect.collidepoint(mouse_pos)
        bg = BUTTON_BG_HOVER if is_hover else BUTTON_BG
        pygame.draw.rect(screen, bg, rect)
        pygame.draw.rect(screen, BLACK, rect, 1)
        label_surf = font.render(label, True, BLACK)
        screen.blit(label_surf, label_surf.get_rect(center=rect.center))

def game(E, Atobe=False, rest=False):
    if Atobe:
        M=conv1(E)
        ir=irarea(E)
        lang=[[-y,x] for [x,y] in rep(E)[0]] 
    else:
        M=E
        ir=[]
        lang=[]
    if M==[]:
        M=randexseg(20,5)
    M=sort(M)
    M=convert_to_int_tuples(M)
    picture(M,lang)
    pygame.init()
    screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
    pygame.display.set_caption("Interval Visualizer")

    font = pygame.font.SysFont(None, 24)
    clock = pygame.time.Clock()

    mode = "main"  

    labels=[]
    for i in range(len(M)):
        labels.append(i+1)
    lang_button = pygame.Rect(LEFT_MARGIN, 10, BUTTON_WIDTH, BUTTON_HEIGHT)
    back_button = pygame.Rect(LEFT_MARGIN + BUTTON_WIDTH + BUTTON_MARGIN, 10, BUTTON_WIDTH, BUTTON_HEIGHT)
    button_rects_main = [("Langlands Data", lang_button), ("Back", back_button)]
    button_rects_lang = [("Back to Main", back_button)]

    while True:
        min_B = min(B for ((A, B), m) in M)
        max_A = max(A for ((A, B), m) in M)
        grid_start = min_B - 1
        grid_end = max_A + 1
        grid_width = grid_end - grid_start + 1
        usable_width = SCREEN_WIDTH - LEFT_MARGIN - RIGHT_MARGIN
        grid_spacing = max(30, usable_width // grid_width)

        if mode == "main":
            diff_boxes, distance_boxes = draw_main_view(screen, M, labels,font, grid_start, grid_end, grid_spacing, button_rects_main)
        else:
            draw_langlands_view(screen, font, button_rects_lang)

        pygame.display.flip()
        clock.tick(60)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                return

            elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                pos = event.pos

                if mode == "main":
                    for (label, rect) in button_rects_main:
                        if rect.collidepoint(pos):
                            if label == "Langlands Data":
                                mode = "langlands"
                    for (rect, i) in diff_boxes:
                        if rect.collidepoint(pos):
                            M,did = swap(M, i,True)
                            if did:
                                labels=sw(labels,i)
                            break
                    for (rect, i) in distance_boxes:
                        if rect.collidepoint(pos):
                            M = ui(M, i)
                            picture(M,lang)
                            break

                elif mode == "langlands":
                    for (label, rect) in button_rects_lang:
                        if rect.collidepoint(pos):
                            if label == "Back to Main":
                                mode = "main"

if __name__ == "__main__":
    main()
